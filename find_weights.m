%% RUN THE WILSON-COWAN MODEL


%% -------------------- parameters --------------------
p = struct;
p.k_E = 0.1603;      p.x0_E = 18;   
p.k_I = 2.1669;      p.x0_I = 0.949;
p.rE  = 0.01820;    p.rI   = 0.05880;
p.tau_E = 2;        p.tau_I = 8;
p.P_E = 16;          p.P_I  = 0.261;   


mode = "IN"; % "GS" for grid search
             % "IN" for inference
T  = 1000;                            % ms
y0 = [0.15; 0.05];                    % initial [E; I]

if mode == "IN"
    p.w_EE = 85;       p.w_IE = 12.5; %p.w_EE / a; %76.2720; %            
    p.w_EI = 85;       p.w_II = 7.5;%p.w_EI / b; %24.7277; %p.w_EI / b;
                      
    opts   = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [t, y] = ode45(@(t, y) wilson_cowan(t, y, p), 0:0.05:T, y0, opts);
    E = y(:,1);   I = y(:,2);
    
    %% ---------- oscillation frequency (second half of the run) ----------
    late = t > T/2;
    x    = E(late) - mean(E(late));
    if max(E(late)) - min(E(late)) > 1e-3
        nCross = sum(x(1:end-1) < 0 & x(2:end) >= 0);
        f_osc  = nCross / ((T/2) / 1000);                 % Hz
        fprintf('Oscillation: %.1f Hz,  E %.3f..%.3f,  I %.3f..%.3f\n', ...
                f_osc, min(E(late)), max(E(late)), min(I(late)), max(I(late)));
    else
        fprintf('No oscillation: steady state E = %.3f, I = %.3f\n', E(end), I(end));
    end
    
    %% ---------- nullclines ----------
    % E-nullcline (dE/dt = 0), solved for I as a function of E
    Eg = linspace(1e-4, 1, 2000)';
    sE = Eg ./ (1 - p.rE*Eg);   okE = sE > 0 & sE < 1;
    I_nullE = (p.w_EE*Eg + p.P_E - p.x0_E - log(sE./(1-sE))/p.k_E) / p.w_EI;
    I_nullE(~okE) = NaN;
    
    % I-nullcline (dI/dt = 0), solved for E as a function of I
    Ig = linspace(1e-4, 1, 2000)';
    sI = Ig ./ (1 - p.rI*Ig);   okI = sI > 0 & sI < 1;
    E_nullI = (p.x0_I + log(sI./(1-sI))/p.k_I + p.w_II*Ig - p.P_I) / p.w_IE;
    E_nullI(~okI) = NaN;
    
    %% ---------- plots ----------
    figure('Name', 'Wilson-Cowan');
    
    subplot(1,2,1);
    plot(t, E, 'b', t, I, 'r', 'LineWidth', 1.2);
    xlabel('t [ms]');  ylabel('activity');  legend('E (PY)', 'I (IN)');
    title('Time course');
    
    subplot(1,2,2);
    plot(Eg, I_nullE, 'b--', E_nullI, Ig, 'r--', 'LineWidth', 1.2);  hold on;
    plot(E(late), I(late), 'k', 'LineWidth', 1.5);
    plot(y0(1), y0(2), 'ko', 'MarkerFaceColor', 'k');
    xlabel('E');  ylabel('I');  xlim([0 1]);  ylim([0 1]);
    legend('dE/dt = 0', 'dI/dt = 0', 'trajectory', 'start');
    title('Phase plane');

elseif mode == "GS"             
    wee_v = linspace(85, 95, 10);        % ranges and resolution are placeholders
    wei_v = linspace(75, 85, 10);
    wie_v = linspace(5 , 15, 5);
    wii_v = linspace(0 , 10, 5);
    
    [C1, C2, C3, C4] = ndgrid(wee_v, wei_v, wie_v, wii_v);
    
    den = p.k_E*C1 - 4.0727;
    ok  = den > 0 & C3 > 0 & (p.k_E*C2 .* (p.k_I*C3) > (p.k_I*C4 + 4.1843) .* den) & C1*1.2709 > (C4*4.2374 + 40.471);
    
    my_file = fopen('combinations.txt', 'w');       % Create or overwrite the file
    combos = [C1(ok) C2(ok) C3(ok) C4(ok)];         % one valid (c1,c2,c3,c4) per row
    fprintf(my_file,'%d of %d grid points satisfy the constraint\n', size(combos,1), numel(ok));
    
    results = nan(size(combos,1), 1);
    
    h = waitbar(0, 'Running iterations...');
    for n = 1:size(combos,1)
        p.w_EE = combos(n,1);  p.w_EI = combos(n,2);  p.w_IE = combos(n,3);  p.w_II = combos(n,4);
        opts   = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
        [t, y] = ode45(@(t, y) wilson_cowan(t, y, p), 0:0.05:T, y0, opts);
        E = y(:,1);   I = y(:,2);

        %% ---------- oscillation frequency (second half of the run) ----------
        late = t > T/2;
        x    = E(late) - mean(E(late));
        if max(E(late)) - min(E(late)) > 1e-3
            nCross = sum(x(1:end-1) < 0 & x(2:end) >= 0);
            f_osc  = nCross / ((T/2) / 1000);                 % Hz
            fprintf(my_file, 'Oscillation: %.1f Hz,  E %.3f..%.3f,  I %.3f..%.3f\n. wEE: %.3f, wIE: %.3f, wEI: %.3f, wII: %.3f\n', ...
                    f_osc, min(E(late)), max(E(late)), min(I(late)), max(I(late)), p.w_EE, p.w_IE, p.w_EI, p.w_II);
        else
            %fprintf('No oscillation: steady state E = %.3f, I = %.3f\n', E(end), I(end));
        end
    waitbar(n/size(combos,1), h, sprintf('Iteration %d / %d', n, size(combos,1)));
    end
    
end
    
%% RUN THE WILSON-COWAN MODEL


%% -------------------- parameters --------------------

p.k_E = 0.162;      p.x0_E = 12;   
p.k_I = 1.699;      p.x0_I = 1.237;
p.rE  = 0.01820;    p.rI   = 0.05880;
p.tau_E = 2;        p.tau_I = 8;
p.P_E = 9;          p.P_I  = -0.72;   

a = 9.7;     % ration w_EE/w_IE
b = 18;     % ratio W_EI/w_II

mode = "GS"; % "GS" for grid search
             % "IN" for inference

if mode == "IN"
    p.w_EE = 91;        p.w_IE = p.w_EE / a;            
    p.w_EI = 82;        p.w_II = p.w_EI / b;
    T  = 1000;                            % ms
    y0 = [0.15; 0.05];                    % initial [E; I]

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
    w_EEs = linspace(0,200,10);
    w_IIs = linspace(0,200,10);
      
    for k = 1:size(w_EEs,2)
        for j = 1:size(w_IIs,2)
            p.w_EE = w_EEs(k);
            p.w_II = w_IIs(j);

            p.w_IE = p.w_EE / a;            
            p.w_II = p.w_EI / b;
    
            opts   = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
            [t, y] = ode45(@(t, y) wilson_cowan(t, y, p), 0:0.05:T, y0, opts);
            E = y(:,1);   I = y(:,2);
    
            %% ---------- oscillation frequency (second half of the run) ----------
            late = t > T/2;
            x    = E(late) - mean(E(late));
            if max(E(late)) - min(E(late)) > 1e-3
                nCross = sum(x(1:end-1) < 0 & x(2:end) >= 0);
                f_osc  = nCross / ((T/2) / 1000);                 % Hz
                fprintf('Oscillation: %.1f Hz,  E %.3f..%.3f,  I %.3f..%.3f\n. wEE: %.3f, wIE: %.3f, wEI: %.3f, wII: %.3f\n', ...
                        f_osc, min(E(late)), max(E(late)), min(I(late)), max(I(late)), p.w_EE, p.w_IE, p.w_EI, p.w_II);
            else
                %fprintf('No oscillation: steady state E = %.3f, I = %.3f\n', E(end), I(end));
            end
    
        end
    end
end
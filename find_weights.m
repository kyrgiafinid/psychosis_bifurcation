%% RUN THE WILSON-COWAN MODEL


%% ---------- parameters ----------
% TEST SET: classic Wilson & Cowan (1972) oscillator, time in ms.
% Known to give a limit cycle -> use it to check that everything works.
p = struct(); 
p.w_EE = 16;   p.w_EI = 12;          % onto E: from E, from I
p.w_IE = 15;   p.w_II = 3;           % onto I: from E, from I
p.P_E  = 1.25; p.P_I  = 0;           % external input
p.k_E  = 1.3;  p.x0_E = 4;           % E sigmoid: slope, threshold
p.k_I  = 2;    p.x0_I = 3.7;         % I sigmoid: slope, threshold
p.rE   = 1;    p.rI   = 1;           % refractory terms
p.tau_E = 8;   p.tau_I = 8;          % time constants (ms)

% YOUR MODEL: uncomment and fill in to replace the test set

%p.P_E  = 0.9;          p.P_I  = -0.07;
p.k_E  = 0.162;         p.x0_E = 35.869;  
p.k_I = 1.699;          p.x0_I = 1.237;
p.rE  = 0.01820*89;     p.rI   = 0.05880;         % r (s) * rref (Hz)
p.tau_E = 2;            p.tau_I = 8;
%p.w_EE = 4;            p.w_II = 7; 
%p.w_EI = 3;            p.w_IE = 12;  
p.w_EE = 4;             p.w_EI = 69.5;   p.w_IE = 15;   p.w_II = 4.8;
p.P_E  = 39.5;          p.P_I  = -1.6;  
p.P_E  = 49;            p.P_I  = -2;
p.w_EE = 4;             p.w_EI = 95;
p.w_IE = 27;            p.w_II = 0.5;
p.tau_sE = 2;           p.tau_sI = 8;
p.P_E = 21.5;           p.P_I = -2.85;
p.w_EE = 68;            p.w_EI = 112;   p.w_IE = 11;   p.w_II = 8.6;


p.k_E  = 0.162;         p.x0_E = 35;           % E sigmoid: slope, threshold
p.k_I  = 1.7;           p.x0_I = 1.237;         % I sigmoid: slope, threshold
p.rE = 1.62;            p.rI = 1.0;
p.P_E = 10;             p.P_I = 0.7;
p.w_EE = 143;           p.w_EI = 113;  p.w_IE = 25;   p.w_II = 7.8;

p.k_E = 0.162;          p.x0_E = 12;        % x0_E lowered -- see below
p.k_I = 1.699;          p.x0_I = 1.237;
p.rE  = 0.01820;    p.rI   = 0.05880;
p.tau_E = 2;            p.tau_I = 8;
p.P_E = 9;              p.P_I  = -0.72;     % your measured I_total


p.w_EE = 91;            p.w_EI = 0;
p.w_IE = 10;            p.w_II = 0;

%p.w_EE = 134;     p.w_EI = 339;
%p.w_IE = 10.8;     p.w_II = 10;

%p.w_IE = 0.5;       p.w_II = 0.5;


T  = 1000;                            % ms
y0 = [0.15; 0.05];                    % initial [E; I]

%{
E0 = y0(1);
I0 = y0(2);
uE = 45.18; uI = 0.01;
w_EEs = linspace(0,300,100);
w_IIs = linspace(0,300,100);


for k = 1:size(w_EEs,2)
    for j = 1:size(w_IIs,2)
        p.w_EE = w_EEs(k);
        p.w_II = w_IIs(j);

        p.w_EI = (p.w_EE*E0 + P_E - uE) / I0;             % from  w_EE*E0 - w_EI*I0 + P_E = uE
        p.w_IE = (uI - P_I + p.w_II*I0) / E0;             % from  w_IE*E0 - w_II*I0 + P_I = uI

        opts   = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
        [t, y] = ode45(@(t, y) wilson_cowan(t, y, p), 0:0.05:T, y0, opts);
        E = y(:,1);   I = y(:,2);

        %% ---------- simulate ----------


        %% ---------- oscillation frequency (second half of the run) ----------
        late = t > T/2;
        x    = E(late) - mean(E(late));
        if max(E(late)) - min(E(late)) > 1e-3
            nCross = sum(x(1:end-1) < 0 & x(2:end) >= 0);
            f_osc  = nCross / ((T/2) / 1000);                 % Hz
            fprintf('Oscillation: %.1f Hz,  E %.3f..%.3f,  I %.3f..%.3f\n', ...
                    f_osc, min(E(late)), max(E(late)), min(I(late)), max(I(late)));
        else
            %fprintf('No oscillation: steady state E = %.3f, I = %.3f\n', E(end), I(end));
        end

    end
end
%}
opts   = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, y] = ode45(@(t, y) wilson_cowan(t, y, p), 0:0.05:T, y0, opts);
E = y(:,1);   I = y(:,2);

%% ---------- simulate ----------


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
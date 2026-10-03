function [nu, t, V] = firing_rate(type, I, T, T_skip)
% FIRING_RATE  Firing rate of a single isolated neuron for a constant input
%
%   nu = firing_rate('IN', 5)          -> rate of IN cell at I = 5 µA/cm²
%   nu = firing_rate('PY', 5)          -> rate of PY cell at I = 5 µA/cm²
%   [nu, t, V] = firing_rate(...)      -> also returns the voltage trace
%
% Inputs:
%   type   : 'PY' or 'IN'
%   I      : injected current (µA/cm²)
%   T      : simulation length (ms)            [default: IN 1000, PY 2000]
%   T_skip : initial transient to ignore (ms)  [default: IN 200,  PY 1000]
%
% Output:
%   nu     : firing rate (Hz)
%   t, V   : time (ms) and membrane potential (mV)

switch upper(type)
    case 'IN'
        if nargin < 3, T = 1000;     end
        if nargin < 4, T_skip = 200; end
        V0 = -65;
        ah = 0.07*exp(-(V0+58)/20);                bh = 1/(exp(-0.1*(V0+28))+1);
        an = -0.01*(V0+34)/(exp(-0.1*(V0+34))-1);  bn = 0.125*exp(-(V0+44)/80);
        y0  = [V0; ah/(ah+bh); an/(an+bn)];
        rhs = @(t, y) rhs_IN(t, y, I);

    case 'PY'
        if nargin < 3, T = 2000;      end
        if nargin < 4, T_skip = 1000; end
        V0 = -70;
        y0  = [V0; 0.5*(1+tanh((V0+2)/21)); 1/(1+exp(-V0/5))];
        rhs = @(t, y) rhs_PY(t, y, I);

    otherwise
        error('type must be ''PY'' or ''IN''');
end

% ---------- integrate ----------
opts   = odeset('RelTol',1e-6, 'AbsTol',1e-8, 'MaxStep',0.1);
[t, y] = ode45(rhs, 0:0.05:T, y0, opts);
V      = y(:,1);

% ---------- count spikes (upward 0 mV crossings after transient) ----------
up   = (V(1:end-1) < 0) & (V(2:end) >= 0);
keep = t(2:end) > T_skip;
nu   = sum(up & keep) / ((T - T_skip)/1000);   % Hz

end


% ================================================================
function dy = rhs_IN(t, y, I)
% Your IN ODE, no synaptic input, plus injected current (Cm = 1)
[dV, dh, dn] = IN_neuron_ode_vec(t, y(1), y(2), y(3), ...
    0, 0, 0, 0, ...      % gAMPA, gNMDAs, gNMDAf, gGABA
    1, 0, 0.05, 0, ...   % X, Y, c, b
    0, 0, 0, 0, ...      % arrivedAMPA, gGABA_source, extSpike, extGaba
    0.05, 100, ...       % dt, bt
    0, 0);               % spike_own, asyncEvent
dy = [dV + I; dh; dn];
end

% ================================================================
function dy = rhs_PY(t, y, I)
% Your PY ODE, no synaptic input, plus injected current (Cm = 1)
[dV, dw, dz] = PY_neuron_ode_vec(t, y(1), y(2), y(3), ...
    0, 0, 0, 0, ...      % gAMPA, gNMDAs, gNMDAf, gGABA
    0, 0, 0, 0, ...      % arrivedAMPA, gGABA_source, extSpike, extGaba
    0.05);               % dt
dy = [dV + I; dw; dz];
end
clear; clc; close all;

% Parameters
Iamp = 1.0;       % external current [uA/cm2]
Ifreq = 0;
p    = [Iamp; Ifreq; 0; 0];

% Initial conditions
V0 = -65;   % chosen resting potential [mV]

alpha_h = 0.07 * exp(-(V0 + 58) / 20);
beta_h  = 1 / (exp(-0.1 * (V0 + 28)) + 1);
h0 = alpha_h / (alpha_h + beta_h);

alpha_n = -0.01 * (V0 + 34) / (exp(-0.1 * (V0 + 34)) - 1);
beta_n  = 0.125 * exp(-(V0 + 44) / 80);
n0 = alpha_n / (alpha_n + beta_n);

x0 = [V0;        % Vm
      h0;     % h
      n0;         % n
      0;         % gAMPA   - no synapse active yet
      0;         % gNMDA_S - no synapse active yet
      0;         % gNMDA_F - no synapse active yet
      0];        % gGABA   - no synapse active yet

% Simulate
T_total      = 1000;   % total simulation time [ms]
T_transient  = 200;    % discard this much at the start [ms]
T_measure    = (T_total - T_transient) / 1000;  % measurement window [s]

I_range_IN = 0:0.01:3;    % sweep range for IN [uA/cm2] -- adjust after first look
spikeThresh_IN = -20;   % mV

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10,'MaxStep', 0.001);
%[t, X]  = ode45(@(t,x) IN_neuron_ode(t,x,p), tspan, x0, options);

rate_IN = zeros(size(I_range_IN));
h = waitbar(0, 'Processing...');
for k = 1:length(I_range_IN)
    p = [I_range_IN(k), 0, 0, 0];
    [t, X] = ode45(@(t,x) IN_neuron_ode(t,x,p), [0 T_total], x0, options);
    rate_IN(k) = count_spike_rate(t, X(:,1), spikeThresh_IN, T_transient, T_measure);
    % Update loading bar
    waitbar(k/length(I_range_IN), h, sprintf('Processing %d/%d', k, length(I_range_IN)));
end

%% Plot
figure('Color','w');
 
plot(I_range_IN, rate_IN, 'r.-', 'LineWidth', 1.2, 'MarkerSize', 10);
xlabel('I_{amp} [\muA/cm^2]'); ylabel('Firing rate [Hz]');
title('IN f-I curve'); grid on;
 
%% Save results for later use (e.g. fitting a gain function F(I))
save('fI_curves.mat', 'I_range_IN', 'rate_IN', 'I_range_IN', 'rate_IN');

%% ------------------------------------------------------------------
function rate = count_spike_rate(t, V, thresh, T_transient, T_measure)
    % Counts upward threshold crossings of V after T_transient, and
    % converts the count to a rate in Hz using the measurement window
    % duration T_measure [s].
    mask = t >= T_transient;
    tSeg = t(mask);
    VSeg = V(mask);
 
    crossings = (VSeg(1:end-1) < thresh) & (VSeg(2:end) >= thresh);
    nSpikes = sum(crossings);
 
    rate = nSpikes / T_measure;
end


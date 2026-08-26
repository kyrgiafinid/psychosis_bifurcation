clear; clc; close all;

% Parameters
Iamp = 50;       % external current [uA/cm2]
Ifreq = 0;
p    = [Iamp; Ifreq; 0; 0];

% Initial conditions
Vm_rest = -70;
w_rest  = 0.5*(1 + tanh((Vm_rest + 2)  /21));
z_rest  = 1  /(1 + exp(-Vm_rest/5));
x0 = [Vm_rest;   % Vm
      w_rest;    % w
      z_rest;    % z
      0;         % gAMPA   - no synapse active yet
      0;         % gNMDA_S - no synapse active yet
      0;         % gNMDA_F - no synapse active yet
      0];        % gGABA   - no synapse active yet

% Simulate
T_total      = 1000;   % total simulation time [ms]
T_transient  = 200;    % discard this much at the start [ms]
T_measure    = (T_total - T_transient) / 1000;  % measurement window [s]

I_range_PY = 0:10:50;    % sweep range for IN [uA/cm2] -- adjust after first look
spikeThresh_PY = -20;   % mV

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10,'MaxStep', 0.001);

rate_PY = zeros(size(I_range_PY));
h = waitbar(0, 'Processing...');
for k = 1:length(I_range_PY)
    p = [I_range_PY(k), 0, 0, 0];
    [t, X] = ode45(@(t,x) PY_neuron_ode(t,x,p), [0 T_total], x0, options);
    rate_PY(k) = count_spike_rate(t, X(:,1), spikeThresh_PY, T_transient, T_measure);
    % Update loading bar
    waitbar(k/length(I_range_PY), h, sprintf('Processing %d/%d', k, length(I_range_PY)));
end

%% Plot
figure('Color','w');
 
plot(I_range_PY, rate_PY, 'r.-', 'LineWidth', 1.2, 'MarkerSize', 10);
xlabel('I_{amp} [\muA/cm^2]'); ylabel('Firing rate [Hz]');
title('IN f-I curve'); grid on;
 
%% Save results for later use (e.g. fitting a gain function F(I))
save('fI_curves.mat', 'I_range_PY', 'rate_PY', 'I_range_PY', 'rate_PY');

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
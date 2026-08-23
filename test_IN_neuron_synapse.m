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
tspan   = [0 1000];     % 1000 ms
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10,'MaxStep', 0.001);
[t, X]  = ode45(@(t,x) IN_neuron_ode(t,x,p), tspan, x0, options);

% Plot
figure;

subplot(3,1,1);
plot(t, X(:,1), 'b-', 'LineWidth', 1.5);
ylabel('V_m [mV]');
title('Pyramidal Neuron - Membrane Potential');
grid on;

subplot(3,1,2);
plot(t, X(:,2), 'r-', 'LineWidth', 1.5);
ylabel('w (recovery)');
grid on;


subplot(3,1,3);
plot(t, X(:,3), 'g-', 'LineWidth', 1.5);
ylabel('z (adaptation)');
xlabel('Time [ms]');
grid on;

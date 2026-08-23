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
tspan   = [0 1000];     % 1000 ms
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10,'MaxStep', 0.001);
[t, X]  = ode45(@(t,x) PY_neuron_ode(t,x,p), tspan, x0, options);

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

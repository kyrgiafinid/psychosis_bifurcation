% --------------- Parameter calculation for Wilson Cowan model ---------------
% Run this script after having run the network simulation script
% Each section calculates a different set of parameters

%% Calculation of refractory period parameters (re ri)
idx_PY      = 1:NPY;
idx_IN      = NPY+1:NPY+NIN;
T_transient = 300; 

T  = rasterT(1:rasterCount); 
ID = rasterID(1:rasterCount);
keep = T > T_transient;
T  = T(keep);
ID = ID(keep);

% Sort by neuron, then by time
[~, order] = sortrows([ID T]);
ID = ID(order);
T  = T(order);
    
% ISIs between consecutive spikes of the SAME neuron
same_cell = diff(ID) == 0;          
ISI       = diff(T);
ISI       = ISI(same_cell);
ISI_cell  = ID([same_cell; false]); 

% Minimum ISI per neuron
N_total = max([idx_PY(:); idx_IN(:)]);
min_ISI = accumarray(ISI_cell, ISI, [N_total 1], @min, Inf);

min_ISI_E = min_ISI(idx_PY);
min_ISI_I = min_ISI(idx_IN);

fired_E = isfinite(min_ISI_E);
fired_I = isfinite(min_ISI_I);

if ~any(fired_E), error('No PY cell fired twice after the transient'); end
if ~any(fired_I), error('No IN cell fired twice after the transient'); end


% Average across neurons, convert to seconds
tau_ref_E = mean(min_ISI_E(fired_E));   % ms
tau_ref_I = mean(min_ISI_I(fired_I));   % ms

r_E = tau_ref_E / 1000;                 % s
r_I = tau_ref_I / 1000;                 % s

fprintf('PY cells used: %d / %d\n', sum(fired_E), numel(idx_PY));
fprintf('IN cells used: %d / %d\n', sum(fired_I), numel(idx_IN));
fprintf('tau_ref_E = %.2f ms  ->  r_E = %.5f s\n', tau_ref_E, r_E);
fprintf('tau_ref_I = %.2f ms  ->  r_I = %.5f s\n', tau_ref_I, r_I);

%% Compute mean extrenal excitation
% before running this run script with zero connections
T_transient = 100;
post = (1:nSteps)' * dt > T_transient;
tt   = (1:nSteps)' * dt;

fprintf('\n                 PY          IN\n');
fprintf('gAMPA   %10.4f  %10.4f   mS/cm^2\n', mean(gAMPA_PY_mean(post)), mean(gAMPA_IN_mean(post)));
fprintf('gNMDA   %10.4f  %10.4f   mS/cm^2\n', mean(gNMDA_PY_mean(post)), mean(gNMDA_IN_mean(post)));
fprintf('gGABA   %10.4f  %10.4f   mS/cm^2\n', mean(gGABA_PY_mean(post)), mean(gGABA_IN_mean(post)));
fprintf('I_exc   %10.3f  %10.3f   uA/cm^2\n', mean(Iexc_PY_mean(post)),  mean(Iexc_IN_mean(post)));
fprintf('I_inh   %10.3f  %10.3f   uA/cm^2\n', mean(Iinh_PY_mean(post)),  mean(Iinh_IN_mean(post)));
fprintf('I_total %10.3f  %10.3f   uA/cm^2\n', ...
        mean(Iexc_PY_mean(post) + Iinh_PY_mean(post)), mean(Iexc_IN_mean(post) + Iinh_IN_mean(post)));

figure('Name', 'Mean synaptic currents');
subplot(2,1,1);
plot(tt, Iexc_PY_mean, 'b', tt, Iinh_PY_mean, 'r', tt, Iexc_PY_mean + Iinh_PY_mean, 'k');
ylabel('\muA/cm^2');  title('PY');  legend('excitatory', 'inhibitory', 'total');
subplot(2,1,2);
plot(tt, Iexc_IN_mean, 'b', tt, Iinh_IN_mean, 'r', tt, Iexc_IN_mean + Iinh_IN_mean, 'k');
ylabel('\muA/cm^2');  xlabel('t [ms]');  title('IN');  legend('excitatory', 'inhibitory', 'total');

%% Compute parameters for sigmoids
% unrelated to simulation, 
I_list_IN = linspace(-1, 3, 50);
I_list_PY = linspace(-1, 60, 150);

nu_IN  = arrayfun(@(I) firing_rate('IN', I), I_list_IN);
nu_PY  = arrayfun(@(I) firing_rate('PY', I), I_list_PY);

[k_I, x0_I, rref_I] = fit_sigmoid(I_list_IN, nu_IN, 3, "IN");    % IN: fit 0–3 uA/cm^2
[k_E, x0_E, rref_E] = fit_sigmoid(I_list_PY, nu_PY, 60, "PY");   % PY: fit up to 60 uA/cm^2
function dY = networkODE(t, Y, I_amp_PY, I_amp_IN, NPY, NIN, prm)
% Continuous dynamics ONLY. No spike-triggered jumps, no /dt_local terms.
% All discrete events (phasic + async release, external Poisson drive,
% X depletion, clamps) are applied in the driver loop after the RK4 step.
%
% prm: struct with fields gNMDA_scale, gGABA_scale
% Y ordering: [PY block (NPY * nStPY states); IN block (NIN * nStIN states)]

nStPY = 7;    % [V, w, gAMPA, gNMDAs, gNMDAf, gGABA, ...]  <-- set to your actual count
nStIN = 11;   % [V, h, n, gAMPA, gNMDAs, gNMDAf, gGABA, X, Y, c, b]

dY = zeros(size(Y));

for i = 1:NPY
    idx = (i-1)*nStPY + (1:nStPY);
    p_i = [I_amp_PY(i), 0, prm.gNMDA_scale, prm.gGABA_scale]; p_i = [0, 0, arrivedAMPA_PY(i), gGABA_source_PY(i), extSpikePY(i), extGabaPY(i), dt, gNMDA_scale, 0, 0, 0];
    dY(idx) = PY_neuron_ode(t, Y(idx), p_i);
end

offset = NPY*nStPY;

for i = 1:NIN
    idx = offset + (i-1)*nStIN + (1:nStIN);
    p_i = [I_amp_IN(i), 0, prm.gNMDA_scale, prm.gGABA_scale];
    dY(idx) = IN_neuron_ode(t, Y(idx), p_i);
end
end
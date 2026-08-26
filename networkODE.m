function dY = networkODE(t, Y, spikeAMPA, spikeGABA, I_amp_PY, I_amp_IN, NPY, NIN)
% spikeAMPA, spikeGABA: (NPY+NIN)x1 logical, ordered [PY block; IN block],
% same ordering as Y. spikeAMPA(i)=1 if neuron i should see an incoming
% excitatory (PY-sourced) spike this step; spikeGABA(i)=1 for inhibitory
% (IN-sourced).

nStatesPerNeuron = 7;
dY = zeros(size(Y));

for i = 1:NPY
    idx = (i-1)*nStatesPerNeuron + (1:nStatesPerNeuron);
    p_i = [I_amp_PY(i), 0, spikeAMPA(i), spikeGABA(i)];
    dY(idx) = PY_neuron_ode(t, Y(idx), p_i);
end

offset = NPY*nStatesPerNeuron;
for i = 1:NIN
    idx = offset + (i-1)*nStatesPerNeuron + (1:nStatesPerNeuron);
    p_i = [I_amp_IN(i), 0, spikeAMPA(NPY+i), spikeGABA(NPY+i)];
    dY(idx) = IN_neuron_ode(t, Y(idx), p_i);
end
end
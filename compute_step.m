function dY = compute_step(t, Y, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, bt_vec, spikedIN, asyncEvent)
    dY = zeros(size(Y));

    % PY: 7 states each, 11-element p (last 3 unused, kept for layout parity)
    for i = 1:NPY
        idx = (i-1)*7 + (1:7);
        p_i = [0, 0, arrivedAMPA_PY(i), gGABA_source_PY(i), extSpikePY(i), extGabaPY(i), dt];
        dY(idx) = PY_neuron_ode(t, Y(idx), p_i);
    end

    % IN: 11 states each, full 11-element p including bt/spike_own_IN/asyncEvent_own
    off = NPY*7;
    for i = 1:NIN
        idx = off + (i-1)*11 + (1:11);
        p_i = [0, 0, arrivedAMPA_IN(i), gGABA_source_IN(i), extSpikeIN(i), extGabaIN(i), dt, bt_vec(i), spikedIN(i), asyncEvent(i)];
        dY(idx) = IN_neuron_ode(t, Y(idx), p_i);
    end
end
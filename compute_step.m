function dY = compute_step(t, Y, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, bt_vec, spikedIN, asyncEvent,PY_idx_all,IN_idx_all)
    persistent p_i_PY p_i_IN dY_template
    if isempty(p_i_PY)
        p_i_PY = zeros(1,7);
        p_i_IN = zeros(1,11);
        dY_template = zeros(size(Y));
    end
    dY = dY_template;


    % PY: 7 states each, 11-element p (last 3 unused, kept for layout parity)
    for i = 1:NPY
        idx = PY_idx_all{i};
        p_i_PY(3) = arrivedAMPA_PY(i);
        p_i_PY(4) = gGABA_source_PY(i);
        p_i_PY(5) = extSpikePY(i);
        p_i_PY(6) = extGabaPY(i);
        p_i_PY(7) = dt;
        dY(idx) = PY_neuron_ode(t, Y(idx), p_i_PY);
    end

    % IN: 11 states each, full 11-element p including bt/spike_own_IN/asyncEvent_own
    for i = 1:NIN
        idx = IN_idx_all{i};
        p_i_IN(3) = arrivedAMPA_IN(i);
        p_i_IN(4) = gGABA_source_IN(i);
        p_i_IN(5) = extSpikeIN(i);
        p_i_IN(6) = extGabaIN(i);
        p_i_IN(7) = dt;
        p_i_IN(8) = bt_vec(i);
        p_i_IN(9) = spikedIN(i);
        p_i_IN(10) = asyncEvent(i);
        dY(idx) = IN_neuron_ode(t, Y(idx), p_i_IN);
    end
end
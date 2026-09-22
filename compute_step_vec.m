function dY = compute_step_vec(t, Y, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, bt_vec, spikedIN, asyncEvent, idxPY, idxIN)
    % idxPY, idxIN: structs with fields .V,.w,.z,.gAMPA,.gNMDAs,.gNMDAf,.gGABA
    % (idxIN also has .X,.Y,.c,.b) -- the strided index vectors, computed
    % ONCE outside the time loop and passed in unchanged every call.

    persistent dY_template
    if isempty(dY_template)
        dY_template = zeros(size(Y));
    end
    dY = dY_template;

    %% ---- PY population: extract state, call vectorized ODE, scatter back ----
    Vm     = Y(idxPY.V);
    w      = Y(idxPY.w);
    z      = Y(idxPY.z);
    gAMPA  = Y(idxPY.gAMPA);
    gNMDAs = Y(idxPY.gNMDAs);
    gNMDAf = Y(idxPY.gNMDAf);
    gGABA  = Y(idxPY.gGABA);

    [dVm, dw, dz, dgAMPA, dgNMDAs, dgNMDAf, dgGABA] = PY_neuron_ode_vec(t, Vm, w, z, gAMPA, gNMDAs, gNMDAf, gGABA, ...
        arrivedAMPA_PY, gGABA_source_PY, extSpikePY, extGabaPY, dt);

    dY(idxPY.V)      = dVm;
    dY(idxPY.w)      = dw;
    dY(idxPY.z)      = dz;
    dY(idxPY.gAMPA)  = dgAMPA;
    dY(idxPY.gNMDAs) = dgNMDAs;
    dY(idxPY.gNMDAf) = dgNMDAf;
    dY(idxPY.gGABA)  = dgGABA;

    %% ---- IN population: same idea, 11 states including GABA terminal ----
    V      = Y(idxIN.V);
    h      = Y(idxIN.h);
    n      = Y(idxIN.n);
    gAMPA_I = Y(idxIN.gAMPA);
    gNMDAs_I = Y(idxIN.gNMDAs);
    gNMDAf_I = Y(idxIN.gNMDAf);
    gGABA_I = Y(idxIN.gGABA);
    X      = Y(idxIN.X);
    Yrel   = Y(idxIN.Y);
    c      = Y(idxIN.c);
    b      = Y(idxIN.b);

    [dV, dh, dn, dgAMPA_I, dgNMDAs_I, dgNMDAf_I, dgGABA_I, dX, dY_rel, dc, db] = IN_neuron_ode_vec(t, V, h, n, ...
        gAMPA_I, gNMDAs_I, gNMDAf_I, gGABA_I, X, Yrel, c, b, ...
        arrivedAMPA_IN, gGABA_source_IN, extSpikeIN, extGabaIN, dt, bt_vec, spikedIN, asyncEvent);

    dY(idxIN.V)      = dV;
    dY(idxIN.h)      = dh;
    dY(idxIN.n)      = dn;
    dY(idxIN.gAMPA)  = dgAMPA_I;
    dY(idxIN.gNMDAs) = dgNMDAs_I;
    dY(idxIN.gNMDAf) = dgNMDAf_I;
    dY(idxIN.gGABA)  = dgGABA_I;
    dY(idxIN.X)      = dX;
    dY(idxIN.Y)      = dY_rel;
    dY(idxIN.c)      = dc;
    dY(idxIN.b)      = db;
end

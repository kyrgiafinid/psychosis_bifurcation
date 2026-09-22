function [dV, dh, dn, dgAMPA, dgNMDAs, dgNMDAf, dgGABA, dX, dY, dc, db] = IN_neuron_ode_vec(t, V, h, n, gAMPA, gNMDAs, gNMDAf, gGABA, X, Y, c, b, arrivedAMPA_IN, gGABA_source_IN, extSpikeIN, extGabaIN, dt, bt_vec, spike_own_IN, asyncEvent_own)
    % Vectorized version of IN_neuron_ode.m -- every input is an NINx1
    % column vector; every output is NINx1. Same equations/constants as
    % the scalar version, elementwise operators throughout.

    persistent gNa gK gL Cm phi EAMPA ENa EK EL EGABA tau_AMPA g_AMPA tau_NMDA_S tau_NMDA_F g_NMDA tau_GABA g_EXT_INPY g_EXT_ININ ...
        tau_D tau_R beta KP kminus kplus c_rest c0 IP U_gaba Util_gaba gama
    if isempty(gNa)
        gNa = 35; gK = 9; gL = 0.1; Cm = 1; phi = 5;
        ENa = 55; EK = -90; EL = -65; EAMPA = 0; EGABA = -75;
        tau_AMPA = 2; g_AMPA = 2e-3;
        tau_NMDA_S = 50; tau_NMDA_F = 2; g_NMDA = 0.1*g_AMPA;
        tau_GABA = 8;
        g_EXT_INPY = 0.003; g_EXT_ININ = 1e-4;
        tau_D = 2; tau_R = 200;
        beta = 5e-3; KP = 0.4; kminus = 0.95e-3; kplus = kminus/0.051;
        c_rest = 0.05; c0 = 20;
        IP = beta*c_rest^2/(c_rest^2+KP^2);
        U_gaba = 0.3; Util_gaba = 0.03; gama = 0.08;
    end

    dt_local = dt;

    %% Gating kinetics (Wang-Buzsaki), elementwise
    alpha_m = -0.1 * (V + 35) ./ (exp(-0.1 * (V + 35)) - 1);
    beta_m  = 4 * exp(-(V + 60) / 18);
    m_inf   = alpha_m ./ (alpha_m + beta_m);

    alpha_h = 0.07 * exp(-(V + 58) / 20);
    beta_h  = 1 ./ (exp(-0.1 * (V + 28)) + 1);

    alpha_n = -0.01 * (V + 34) ./ (exp(-0.1 * (V + 34)) - 1);
    beta_n  = 0.125 * exp(-(V + 44) / 80);

    %% Ionic currents
    INa = gNa * m_inf.^3 .* h .* (V - ENa);
    IK  = gK  * n.^4       .* (V - EK);
    IL  = gL              .* (V - EL);
    Iion = INa + IK + IL;

    %% NMDA voltage dependence
    gNMDA_total = (gNMDAs - gNMDAf) ./ (1 + 0.264 * exp(-0.06 * V));

    %% Synaptic currents
    Iampa = gAMPA       .* (V - EAMPA);
    Inmda = gNMDA_total .* (V - EAMPA);
    Igaba = gGABA       .* (V - EGABA);
    Isyn = Iampa + Inmda + Igaba;

    %% Membrane potential ODE
    dV = (-Iion - Isyn) / Cm;

    %% Gating ODEs
    dh = phi * (alpha_h .* (1 - h) - beta_h .* h);
    dn = phi * (alpha_n .* (1 - n) - beta_n .* n);

    %% Synaptic conductance ODEs
    dgAMPA  = -gAMPA/tau_AMPA    + (g_AMPA/dt_local)*arrivedAMPA_IN      + (g_EXT_INPY/dt_local)*extSpikeIN;
    dgNMDAs = -gNMDAs/tau_NMDA_S + (g_NMDA/dt_local)*arrivedAMPA_IN      + (0.1*g_EXT_INPY/dt_local)*extSpikeIN;
    dgNMDAf = -gNMDAf/tau_NMDA_F + (g_NMDA/dt_local)*arrivedAMPA_IN      + (0.1*g_EXT_INPY/dt_local)*extSpikeIN;
    dgGABA  = -gGABA/tau_GABA    + gGABA_source_IN                       + (g_EXT_ININ/dt_local)*extGabaIN;

    %% GABA release terminal, fully continuous (jumps folded in via /dt trick)
    Z = 1 - X - Y;
    dX = Z/tau_R - (U_gaba/dt_local)*X.*spike_own_IN - (Util_gaba/dt_local)*X.*asyncEvent_own;
    dY = -Y/tau_D + (U_gaba/dt_local)*X.*spike_own_IN + (Util_gaba/dt_local)*X.*asyncEvent_own;
    pump = beta*c.^2 ./ (c.^2 + KP^2);
    buffer = kplus*(bt_vec-b).*c - kminus*b;

    c_safe = max(c, 1e-6);
    dc_spike = gama * log(c0./c_safe);

    dc = -pump - buffer + IP + (dc_spike/dt_local).*spike_own_IN;
    db = buffer;
end

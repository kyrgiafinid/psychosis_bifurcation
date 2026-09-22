function [dVm, dw, dz, dgAMPA, dgNMDAs, dgNMDAf, dgGABA] = PY_neuron_ode_vec(t, Vm, w, z, gAMPA, gNMDAs, gNMDAf, gGABA, arrivedAMPA_PY, gGABA_source_PY, extSpikePY, extGabaPY, dt)
    % Vectorized version of PY_neuron_ode.m -- every input is an NPYx1
    % column vector (one entry per PY neuron); every output is NPYx1 too.
    % Same equations, same constants, as the scalar version -- only the
    % arithmetic operators changed to element-wise (.*  ./  .^) so they
    % apply to the whole population at once instead of looping.

    persistent GNa GK GL GA Cm EAMPA EGABA tau_AMPA tau_NMDA_S tau_NMDA_F g_NMDA g_AMPA tau_GABA g_EXT_PYPY g_EXT_PYIN
    if isempty(GNa)
        GNa = 10; GK = 10; GL = 1.3; GA = 3; Cm = 1;
        EAMPA = 0; EGABA = -75;
        tau_AMPA = 2; g_AMPA = 7.5e-3;
        tau_NMDA_S = 100; tau_NMDA_F = 2; g_NMDA = 0.4*g_AMPA;
        tau_GABA = 8;
        g_EXT_PYPY = 0.25; g_EXT_PYIN = 0.025;
    end

    dt_local = dt;   % scalar, same for every neuron

    %% Gating variables (elementwise)
    m_inf = 0.5 * (1 + tanh((Vm + 1.2) / 23));
    w_inf = 0.5 * (1 + tanh((Vm + 2)   / 21));
    z_inf = 1   ./ (1 + exp(-Vm/5));

    %% Ionic current
    Iion = GNa * m_inf .* (Vm - 50) ...
         + GK  * w      .* (Vm + 100) ...
         + GL           .* (Vm + 70) ...
         + GA  * z      .* (Vm + 100);

    %% NMDA voltage dependence (Mg2+ block)
    gNMDA_total = (gNMDAs - gNMDAf) ./ (1 + 0.264*exp(-0.06*Vm));

    %% Synaptic currents
    Iampa = gAMPA       .* (Vm - EAMPA);
    Inmda = gNMDA_total .* (Vm - EAMPA);
    Igaba = gGABA       .* (Vm - EGABA);
    Isyn = Iampa + Inmda + Igaba;

    %% Membrane potential ODE
    dVm = (-Iion - Isyn) / Cm;

    %% Gating ODEs
    dw = 0.15 * (w_inf - w) .* cosh((Vm + 2) / 42);
    dz = 0.005 * (z_inf - z);

    %% Synaptic conductance ODEs
    dgAMPA  = -gAMPA/tau_AMPA    + (g_AMPA/dt_local)*arrivedAMPA_PY      + (g_EXT_PYPY/dt_local)*extSpikePY;
    dgNMDAs = -gNMDAs/tau_NMDA_S + (g_NMDA/dt_local)*arrivedAMPA_PY      + (0.4*g_EXT_PYPY/dt_local)*extSpikePY;
    dgNMDAf = -gNMDAf/tau_NMDA_F + (g_NMDA/dt_local)*arrivedAMPA_PY      + (0.4*g_EXT_PYPY/dt_local)*extSpikePY;
    dgGABA  = -gGABA/tau_GABA    + gGABA_source_PY                       + (g_EXT_PYIN/dt_local)*extGabaPY;
end

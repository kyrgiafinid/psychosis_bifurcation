function f = PY_neuron_ode(t, x, p)
    % Pyramidal neuron with synaptic currents
    % t unused - autonomous system
    %
    % States:
    %   x(1) = Vm      - membrane potential [mV]
    %   x(2) = w       - K recovery variable
    %   x(3) = z       - adaptation variable
    %   x(4) = gAMPA   - AMPA conductance
    %   x(5) = gNMDA_S - slow NMDA conductance
    %   x(6) = gNMDA_F - fast NMDA conductance
    %   x(7) = gGABA   - GABA conductance
    %
    % Parameters:
    %   p(1) = Iext             - external current
    %   p(2) = I_freq           - external exitation frequency
    %   p(3) = spike_AMPA       - presynaptic AMPA spike train (0 or 1)
    %   p(4) = spike_GABA       - presynaptic GABA spike train (0 or 1)
    %   p(5) = spike_AMPA_ext   - binary external exitatory arival
    %   p(6) = spike_GABA_ext   - binary external inhibitory arival
    %   p(7) = dt_local
    
    %% Unpack states

    % Time varying input
    I_amp  = p(1);   % amplitude
    I_freq = p(2);   % frequency [Hz]
    
    % Sinusoidal drive
    Iext = I_amp * (1 + sin(2*pi*I_freq*t/1000));

    Vm     = x(1);
    w      = x(2);
    z      = x(3);
    gAMPA  = x(4);
    gNMDAs = x(5);
    gNMDAf = x(6);
    gGABA  = x(7);
    
    %% Unpack parameters
    spike_AMPA      = p(3);   % 1 if presynaptic excitatory spike, else 0
    gGABA_source_PY      = p(4);   % 1 if presynaptic inhibitory spike, else 0
    spike_AMPA_ext  = p(5);   % external excitatory arrival (binary)
    spike_GABA_ext  = p(6);   % external inhibitory arrival (binary)
    dt_local        = p(7);
    
    %% Membrane parameters
    GNa = 10;
    GK  = 10;
    GL  = 1.3;
    GA  = 3;
    Cm  = 1;
    
    %% Reversal potentials
    EAMPA = 0;    % excitatory reversal
    EGABA = -75;  % inhibitory reversal 
    
    %% Synaptic conductance parameters 
    % AMPA
    tau_AMPA  = 2;      % decay time [ms]
    g_AMPA    = 7.5e-3; % peak conductance [mS/cm2]
    
    % NMDA (two components: slow and fast)
    tau_NMDA_S = 100;   % slow decay [ms]
    tau_NMDA_F = 2;     % fast decay [ms]
    g_NMDA     = 0.4 * g_AMPA;  % 40% of AMPA
    
    % GABA
    tau_GABA  = 8;      % decay time [ms]
    g_GABA    = 0.8e-3; % peak conductance [mS/cm2]

    g_EXT_PYPY = 0.25;
    g_EXT_PYIN = 0.025;
    
    %% Gating variables
    m_inf = 0.5 * (1 + tanh((Vm + 1.2) / 23));
    w_inf = 0.5 * (1 + tanh((Vm + 2)   / 21));
    z_inf = 1   / (1 + exp(-Vm/5));
    
    %% Ionic current
    Iion = GNa * m_inf * (Vm - 50) ...      % Na: fast m_inf drives upstroke
     + GK  * w     * (Vm + 100)  ...        % K:  slow w drives repolarization
     + GL          * (Vm + 70)  ...
     + GA  * z     * (Vm + 100);
    %% NMDA voltage dependence (Mg2+ block, from paper eq 10)
    % gNMDA_total = (gNMDAs + gNMDAf) / (1 + 0.264*exp(-0.06*Vm))
    gNMDA_total = (gNMDAs - gNMDAf) / (1 + 0.264*exp(-0.06*Vm));
    
    %% Synaptic currents
    Iampa = gAMPA      * (Vm - EAMPA);   % excitatory
    Inmda = gNMDA_total* (Vm - EAMPA);   % excitatory
    Igaba = gGABA      * (Vm - EGABA);   % inhibitory
    
    % Total synaptic current
    % AMPA and NMDA are excitatory (depolarizing, negative sign removes them)
    % GABA is inhibitory
    Isyn = Iampa + Inmda + Igaba;

    %fprintf('PY Isyn = %.5f', Isyn)
    
    %% Membrane potential ODE
    dVm = (-Iion - Isyn + Iext) / Cm;
    
    %% Gating ODEs
    dw = 0.15 * (w_inf - w) * cosh((Vm - 2) / 42);
    dz = 0.005 * (z_inf - z);
    
    scale_AMPA = 5;
    %% Synaptic conductance ODEs (eq 8-9 from paper)
    % Each conductance decays exponentially and jumps at each spike
    dgAMPA  = -gAMPA/tau_AMPA + (g_AMPA/dt_local)*spike_AMPA * scale_AMPA + (g_EXT_PYPY/dt_local)*spike_AMPA_ext;
    dgNMDAs = -gNMDAs/tau_NMDA_S + (g_NMDA/dt_local)*spike_AMPA* scale_AMPA + (0.4*g_EXT_PYPY/dt_local)*spike_AMPA_ext;
    dgNMDAf = -gNMDAf/tau_NMDA_F + (g_NMDA/dt_local)*spike_AMPA* scale_AMPA + (0.4*g_EXT_PYPY/dt_local)*spike_AMPA_ext;
    dgGABA  = 0; %-gGABA/tau_GABA + gGABA_source_PY + (g_EXT_PYIN/dt_local)*spike_GABA_ext;

    %if gGABA_source_PY ~= 0
    %    gGABA_source_PY
    %end
    if spike_AMPA ~= 0
        spike_AMPA;
        dt_local;
    end
    
    %% Pack output
    f = [dVm; dw; dz; dgAMPA; dgNMDAs; dgNMDAf; dgGABA];
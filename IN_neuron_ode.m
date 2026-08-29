function f = IN_neuron_ode(t, x, p)
    % Interneuronneuron with synaptic currents
    % t unused - autonomous system
    %
    % States:
    %   x(1) = V       - membrane potential [mV]
    %   x(2) = h       - Na inactivation variable
    %   x(3) = n       - K activation variable
    %   x(4) = gAMPA   - AMPA conductance
    %   x(5) = gNMDA_S - slow NMDA conductance
    %   x(6) = gNMDA_F - fast NMDA conductance
    %   x(7) = gGABA   - GABA conductance
    %
    % Parameters:
    %   p(1) = Iext       - external current
    %   p(2) = I_freq     - frequency [Hz] of external drive
    %   p(3) = spike_AMPA - presynaptic PY spike train (0 or 1), drives AMPA/NMDA
    %   p(4) = spike_GABA - presynaptic IN spike train (0 or 1), drives GABA
    %   p(5) = spike_AMPA_ext   - binary external exitatory arival
    %   p(6) = spike_GABA_ext   - binary external inhibitory arival
    %   p(7) = dt_local
    
    %% Unpack states

    % Time varying input
    I_amp  = p(1);   % amplitude
    I_freq = p(2);   % frequency [Hz]
    
    % Sinusoidal drive
    Iext = I_amp * (1 + sin(2*pi*I_freq*t/1000));

    V     = x(1);
    h      = x(2);
    n      = x(3);
    gAMPA  = x(4);
    gNMDAs = x(5);
    gNMDAf = x(6);
    gGABA  = x(7);
    
    %% Unpack parameters
    spike_AMPA = p(3);   % 1 if presynaptic excitatory spike, else 0
    spike_GABA = p(4);   % 1 if presynaptic inhibitory spike, else 0
    spike_AMPA_ext = p(5);
    spike_GABA_ext = p(6);
    dt_local       = p(7);
    
    %% Membrane parameters
    gNa = 35;
    gK  = 9;
    gL  = 0.1;
    Cm  = 1;
    phi = 5;
    
    %% Reversal potentials
    ENa   = 55;    % Na reversal
    EK    = -90;   % K reversal
    EL    = -65;   % leak reversal
    EAMPA = 0;     % excitatory reversal
    EGABA = -75;   % inhibitory reversal
    
    %% Synaptic conductance parameters (PY -> IN and IN -> IN, from paper)
    % AMPA (from PY cells)
    tau_AMPA  = 2;        % decay time [ms]
    g_AMPA    = 2e-3;     % peak conductance [mS/cm2]  (g_AMPA^{IN<-PY})
 
    % NMDA (from PY cells; faster slow-component decay and smaller peak than PY-PY, eq 11)
    tau_NMDA_S = 50;      % slow decay [ms] = 0.5 * tau_NMDA_S^{PY<-PY}
    tau_NMDA_F = 2;       % fast decay [ms]
    g_NMDA     = 0.1 * g_AMPA;  % 10% of AMPA (g_NMDA^{IN<-PY})
 
    % GABA (from other IN cells)
    tau_GABA  = 8;        % decay time [ms]
    g_GABA    = 5e-4;     % peak conductance [mS/cm2]  (g_GABA^{IN<-IN})

    g_EXT_INPY = 0.003;
    g_EXT_ININ = 1e-4;
    
    %% Gating kinetics (Wang-Buzsaki)
    alpha_m = -0.1 * (V + 35) ./ (exp(-0.1 * (V + 35)) - 1);
    beta_m  = 4 * exp(-(V + 60) / 18);
    m_inf   = alpha_m / (alpha_m + beta_m);
 
    alpha_h = 0.07 * exp(-(V + 58) / 20);
    beta_h  = 1 / (exp(-0.1 * (V + 28)) + 1);
 
    alpha_n = -0.01 * (V + 34) ./ (exp(-0.1 * (V + 34)) - 1);
    beta_n  = 0.125 * exp(-(V + 44) / 80);
 
    %% Ionic currents
    INa = gNa * m_inf^3 * h * (V - ENa);
    IK  = gK  * n^4      * (V - EK);
    IL  = gL             * (V - EL);
    Iion = INa + IK + IL;
    
    %% NMDA voltage dependence (Mg2+ block, from paper eq 10)
    gNMDA_total = (gNMDAs - gNMDAf) / (1 + 0.264 * exp(-0.06 * V));
 
    %% Synaptic currents
    Iampa = gAMPA       * (V - EAMPA);   % excitatory (from PY)
    Inmda = gNMDA_total * (V - EAMPA);   % excitatory (from PY)
    Igaba = gGABA       * (V - EGABA);   % inhibitory (from other IN)
 
    Isyn = Iampa + Inmda + Igaba;
 
    %% Membrane potential ODE
    dV = (-Iion - Isyn + Iext) / Cm;
 
    %% Gating ODEs
    dh = phi * (alpha_h * (1 - h) - beta_h * h);
    dn = phi * (alpha_n * (1 - n) - beta_n * n);
 
    %% Synaptic conductance ODEs (eq 8-9 from paper)
    dgAMPA  = -gAMPA/tau_AMPA + (g_AMPA/dt_local)*spike_AMPA + (g_EXT_INPY/dt_local)*spike_AMPA_ext;
    dgNMDAs = -gNMDAs/tau_NMDA_S + (g_NMDA/dt_local)*spike_AMPA + (0.1*g_EXT_INPY/dt_local)*spike_AMPA_ext;
    dgNMDAf = -gNMDAf/tau_NMDA_F + (g_NMDA/dt_local)*spike_AMPA + (0.1*g_EXT_INPY/dt_local)*spike_AMPA_ext;
    dgGABA  = -gGABA/tau_GABA + (g_GABA/dt_local)*spike_GABA + (g_EXT_ININ/dt_local)*spike_GABA_ext;
    
    if spike_GABA ~= 0
        spike_GABA
    end
    if spike_AMPA ~= 0
        spike_AMPA
    end
 
    %% Pack output
    f = [dV; dh; dn; dgAMPA; dgNMDAs; dgNMDAf; dgGABA];
 
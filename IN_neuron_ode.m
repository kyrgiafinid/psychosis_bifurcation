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
    %   x(8)=X  x(9)=Y  x(10)=c  x(11)=b     <- GABA release terminal state
    %
    % Parameters:
    %   p(1) = Iext       - external current
    %   p(2) = I_freq     - frequency [Hz] of external drive
    %   p(3) = spike_AMPA - presynaptic PY spike train (0 or 1), drives AMPA/NMDA
    %   p(4) = spike_GABA - presynaptic IN spike train (0 or 1), drives GABA
    %   p(5) = spike_AMPA_ext   - binary external exitatory arival
    %   p(6) = spike_GABA_ext   - binary external inhibitory arival
    %   p(7) = dt_local
    %   p(8) = gNMDA_scale
    %   p(9) = bt (PV concentration for THIS neuron's terminal)
    %   p(10)= spike_own_IN (this neuron's own spike, triggers phasic release)
    %   p(11)= asyncEvent_own (stochastic async-release flag, decided before RK4)
    
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
    X = x(8); Y = x(9); c = x(10); b = x(11);
    
    %% Unpack parameters
    spike_AMPA = p(3);   % 1 if presynaptic excitatory spike, else 0
    gGABA_source_IN = p(4);   % 1 if presynaptic inhibitory spike, else 0
    spike_AMPA_ext = p(5);
    spike_GABA_ext = p(6);
    dt_local       = p(7);
    gNMDA_scale     = p(8);
    bt              = p(9);
    spike_own_IN    = p(10);   % this neuron's OWN spike this step (phasic trigger)
    asyncEvent_own  = p(11);   % stochastic async-release flag, decided before RK4, frozen across substeps
    
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
    g_AMPA    = 2e-3;     % peak conductance [mS/cm2] 
 
    % NMDA (from PY cells; faster slow-component decay and smaller peak than PY-PY, eq 11)
    tau_NMDA_S = 50;      % slow decay [ms] = 0.5 * tau_NMDA_S^{PY<-PY}
    tau_NMDA_F = 2;       % fast decay [ms]
    g_NMDA     = gNMDA_scale * 0.1 * g_AMPA;  % 10% of AMPA (g_NMDA^{IN<-PY})
 
    % GABA (from other IN cells)
    tau_GABA  = 8;        % decay time [ms]

    g_EXT_INPY = gNMDA_scale * 0.003;
    g_EXT_ININ = 1e-4;

    % ---- GABA release terminal constants (Tsodyks-Markram + calcium/PV) ----
    tau_D=2; tau_R=200;
    alpha_ca=5e-3; KP=0.4; kminus=0.95e-3; kplus=kminus/0.051;
    c_rest=0.05;
    IP = alpha_ca*c_rest^2/(c_rest^2+KP^2);

    %tau_D=2; tau_R=200; kminus = 0.95; kplus=kminus/0.051; 
    %alpha_ca=5e-3; KP=0.4; c_rest=0.05;
    %IP = alpha_ca*c_rest^2/(c_rest^2+KP^2);
    
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

    %fprintf('IN Isyn = %.5f', Isyn)
 
    %% Membrane potential ODE
    dV = (-Iion - Isyn + Iext) / Cm;
 
    %% Gating ODEs
    dh = phi * (alpha_h * (1 - h) - beta_h * h);
    dn = phi * (alpha_n * (1 - n) - beta_n * n);
    
    %% Synaptic conductance ODEs (eq 8-9 from paper)
    dgAMPA  = -gAMPA/tau_AMPA + (g_AMPA/dt_local)*spike_AMPA + (g_EXT_INPY/dt_local)*spike_AMPA_ext;
    dgNMDAs = -gNMDAs/tau_NMDA_S + (g_NMDA/dt_local)*spike_AMPA + (0.1*g_EXT_INPY/dt_local)*spike_AMPA_ext;
    dgNMDAf = -gNMDAf/tau_NMDA_F + (g_NMDA/dt_local)*spike_AMPA + (0.1*g_EXT_INPY/dt_local)*spike_AMPA_ext;
    dgGABA  = -gGABA/tau_GABA + gGABA_source_IN + (g_EXT_ININ/dt_local)*spike_GABA_ext;
    
    U_gaba=0.3; Util_gaba=0.03; beta_ca=0.08;
    
    % ---- GABA release terminal, now FULLY continuous (jumps folded in via /dt trick) ----
    Z = 1 - X - Y;
    dX = Z/tau_R - (U_gaba/dt_local)*X*spike_own_IN - (Util_gaba/dt_local)*X*asyncEvent_own;
    dY = -Y/tau_D + (U_gaba/dt_local)*X*spike_own_IN + (Util_gaba/dt_local)*X*asyncEvent_own;
    pump = alpha_ca*c^2/(c^2+KP^2);
    buffer = kplus*(bt-b)*c - kminus*b;

    % dc spike term, paper eq 17
    c_safe = max(c, 1e-6);                     % guard the log at c -> 0
    dc_spike = beta_ca * log(c_rest/c_safe);
    dc_spike = max(dc_spike, 0);    

    dc = -pump - buffer + IP + (dc_spike/dt_local)*spike_own_IN;
    db = buffer;


    %% Pack output
    f = [dV; dh; dn; dgAMPA; dgNMDAs; dgNMDAf; dgGABA; dX; dY; dc; db];
end
 
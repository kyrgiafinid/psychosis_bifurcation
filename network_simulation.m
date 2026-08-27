clear;
clc;
close all;

%% ============================================================
%  NETWORK PARAMETERS
%  ============================================================

Nx = 30;
Ny = 30;

N = Nx * Ny;
gGABA_scale = 8;
gNMDA_scale = 1;
integration_method = "R4K"; % R4K, Euler
gaba_simulation = "simple"; % simple, ode

rng(1);               % Fixed seed for random number generation

% Neuron types
% 1 = PY
% 2 = IN

type = ones(N,1);
type(5:5:N) = 2;

PY = find(type == 1);
IN = find(type == 2);

NPY = length(PY);
NIN = length(IN);

fprintf('Total neurons: %d\n',N);
fprintf('PY neurons:     %d\n',NPY);
fprintf('IN neurons:     %d\n',NIN);


%% ============================================================
%  LATTICE
%  ============================================================

[xGrid,yGrid] = meshgrid(1:Nx,1:Ny);

x = xGrid(:);
y = yGrid(:);


%% Plot lattice

figure;

scatter(x(PY),y(PY),50,'b','filled');
hold on;

scatter(x(IN),y(IN),70,'r','filled');

xlabel('x');
ylabel('y');

title('8 x 8 neuronal lattice');

legend('PY','IN');

axis equal;
xlim([0.5 Nx+0.5]);
ylim([0.5 Ny+0.5]);

grid on;


%% ============================================================
%  CONNECTIVITY PARAMETERS
%  ============================================================

% Connection probabilities

p_PY_PY = 0.4;
p_PY_IN = 0.3;
p_IN_PY = 0.6;
p_IN_IN = 0.7;

% Spatial footprints

LF_PY_PY = 10;
LF_PY_IN = 10;
LF_IN_PY = 20;

% For IN -> IN the paper specifies probability,
% but does not give a separate footprint.
LF_IN_IN = 10;


%% ============================================================
%  CREATE CONNECTIVITY MATRIX
%  ============================================================

% C(post,pre) = 1 means:
%
%       pre  --->  post

C = zeros(N,N);


for pre = 1:N
    for post = 1:N
        % No self-connections
        if pre == post
            continue;
        end

        % Determine connection type
        if type(pre) == 1 && type(post) == 1
            % PY -> PY
            p = p_PY_PY;
            LF = LF_PY_PY;
        elseif type(pre) == 1 && type(post) == 2
            % PY -> IN
            p = p_PY_IN;
            LF = LF_PY_IN;
        elseif type(pre) == 2 && type(post) == 1
            % IN -> PY
            p = p_IN_PY;
            LF = LF_IN_PY;
        else
            % IN -> IN
            p = p_IN_IN;
            LF = LF_IN_IN;
        end


        % Check spatial footprint
        dx = abs(x(pre) - x(post));
        dy = abs(y(pre) - y(post));

        insideFootprint = (dx <= LF) && (dy <= LF);

        % Randomly create synapse
        if insideFootprint
            if rand < p
                C(post,pre) = 1;
            end
        end
    end
end


%% ============================================================
%  DISPLAY CONNECTIVITY
%  ============================================================

figure;

imagesc(C);

xlabel('Presynaptic neuron');
ylabel('Postsynaptic neuron');

title('Connectivity matrix');

colorbar;


%% ============================================================
%  NUMBER OF CONNECTIONS
%  ============================================================

fprintf('\nConnections:\n');

fprintf('PY -> PY: %d\n', sum(C(PY,PY),'all'));
fprintf('PY -> IN: %d\n', sum(C(IN,PY),'all'));
fprintf('IN -> PY: %d\n', sum(C(PY,IN),'all'));
fprintf('IN -> IN: %d\n',  sum(C(IN,IN),'all'));


%% ============================================================
%  INITIAL STATE VECTOR
%  ============================================================

% Parameters for Tsodyks-Markram + calcium/PV GABA release parameters (paper eqs. 13-18)
U_gaba    = 0.3;
Util_gaba = 0.03;
tau_D     = 2;
tau_R     = 200;
lambda_t  = 0.03e-3;
Ka        = 0.2;
alpha_ca  = 5e-3;
KP        = 0.4;
beta_ca   = 0.08;
kplus     = 0.95e-3;
kminus    = kplus/0.051;
bt        = 100;
c_rest    = 0.05;
IP        = alpha_ca * c_rest^2 / (c_rest^2 + KP^2);

T_X = ones(NIN,1);
T_Y = zeros(NIN,1);
T_c = c_rest*ones(NIN,1);
T_b = zeros(NIN,1);

% Initial conditions for IN neurons
V0 = -65;   % chosen resting potential [mV]

alpha_h = 0.07 * exp(-(V0 + 58) / 20);
beta_h  = 1 / (exp(-0.1 * (V0 + 28)) + 1);
h0 = alpha_h / (alpha_h + beta_h);

alpha_n = -0.01 * (V0 + 34) / (exp(-0.1 * (V0 + 34)) - 1);
beta_n  = 0.125 * exp(-(V0 + 44) / 80);
n0 = alpha_n / (alpha_n + beta_n);

Y_IN = [V0;         % Vm
      h0;           % h
      n0;           % n
      0;            % gAMPA   - no synapse active yet
      0;            % gNMDA_S - no synapse active yet
      0;            % gNMDA_F - no synapse active yet
      0];           % gGABA   - no synapse active yet


% Initial conditions for IN neurons
Vm_rest = -70;
w_rest  = 0.5*(1 + tanh((Vm_rest + 2)  /21));
z_rest  = 1  /(1 + exp(-Vm_rest/5));
Y_PY = [Vm_rest;    % Vm
      w_rest;       % w
      z_rest;       % z
      0;            % gAMPA   - no synapse active yet
      0;            % gNMDA_S - no synapse active yet
      0;            % gNMDA_F - no synapse active yet
      0];           % gGABA   - no synapse active yet

% ---- reorder connectivity to match Y's [PY block; IN block] layout ----
C_PYPY = C(PY,PY);   % PY - PY connections
C_INPY = C(PY,IN);   % PY - IN connections
C_PYIN = C(IN,PY);   % IN - PY connections
C_ININ = C(IN,IN);   % IN - IN connections

% Uncomment to add external input
% I_amp_PY = 80 + 5*randn(NPY,1);   
% I_amp_IN = 80 + 5*randn(NIN,1);

I_amp_PY = zeros(NPY,1);
I_amp_IN = zeros(NIN,1);

dt = 0.05;          
T_total = 500;
nSteps = round(T_total/dt);

Y = [repmat(Y_PY, NPY, 1); repmat(Y_IN, NIN, 1)];
Vidx_PY = 1:7:7*NPY;
Vidx_IN = 7*NPY+1:7:7*NPY+7*NIN;

rasterT = []; rasterID = [];
Vprev = Y([Vidx_PY, Vidx_IN]);   % for threshold-crossing detection

spikeAMPA = false(NPY+NIN,1);
spikeGABA = false(NPY+NIN,1);

LFP = zeros(nSteps,1);

h = waitbar(0, 'Processing...');

% External drives
vSTIM_PY = 250;
p_ext_PY = vSTIM_PY * dt / 1000;

vSTIM_IN = 500;                    
p_ext_IN = vSTIM_IN * dt / 1000;

g_EXT_PYPY = 0.25;    % paper's value, PY's external excitatory conductance
g_EXT_INPY = 0.003;   % paper's value, IN's external excitatory conductance
g_EXT_PYIN = 0.025;   % external inhibitory conductance onto PY (g~EXT^{PY<-IN})
g_EXT_ININ = 1e-4;    % external inhibitory conductance onto IN (g~EXT^{IN<-IN}

gAMPA_idx_PY = 4:7:7*NPY;
gAMPA_idx_IN = 7*NPY+4:7:7*NPY+7*NIN;
gGABA_idx_PY = 7:7:7*NPY;
gGABA_idx_IN = 7*NPY+7:7:7*NPY+7*NIN;
gNMDAs_idx_PY = 5:7:7*NPY;  gNMDAf_idx_PY = 6:7:7*NPY;
gNMDAs_idx_IN = 7*NPY+5:7:7*NPY+7*NIN;  gNMDAf_idx_IN = 7*NPY+6:7:7*NPY+7*NIN;

% testing plots 
zTrace = zeros(nSteps,1);
gGABATrace = zeros(nSteps,1);

for k = 1:nSteps
    t = k*dt;
    extSpikePY = rand(NPY,1) < p_ext_PY;
    extSpikeIN = rand(NIN,1) < p_ext_IN;
    extGabaPY = rand(NPY,1) < p_ext_PY;   % same rate as extSpikePY, independent draw
    extGabaIN = rand(NIN,1) < p_ext_IN;   % same rate as extSpikeIN, independent draw

    % ---- 1. decay-only integration ----
    if integration_method == "Euler"
        dY = networkODE(t, Y, false(NPY+NIN,1), false(NPY+NIN,1), I_amp_PY, I_amp_IN, NPY, NIN);
        Y = Y + dt*dY;
    elseif integration_method == "R4K"
        if gaba_simulation == "ode"
            Yrel = T_Y;   % freeze this step's release fraction for the RK4 substeps
            k1 = computeDY(t,        Y,           NPY, NIN, C_INPY, C_ININ, gGABA_scale*0.8e-3, gGABA_scale*5e-4, Yrel);
            k2 = computeDY(t + dt/2, Y + dt/2*k1, NPY, NIN, C_INPY, C_ININ, gGABA_scale*0.8e-3, gGABA_scale*5e-4, Yrel);
            k3 = computeDY(t + dt/2, Y + dt/2*k2, NPY, NIN, C_INPY, C_ININ, gGABA_scale*0.8e-3, gGABA_scale*5e-4, Yrel);
            k4 = computeDY(t + dt,   Y + dt*k3,   NPY, NIN, C_INPY, C_ININ, gGABA_scale*0.8e-3, gGABA_scale*5e-4, Yrel);
            Y = Y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        elseif gaba_simulation == "simple"
            k1 = computeDY_simple(t,        Y,          NPY, NIN);
            k2 = computeDY_simple(t + dt/2, Y + dt/2*k1, NPY, NIN);
            k3 = computeDY_simple(t + dt/2, Y + dt/2*k2, NPY, NIN);
            k4 = computeDY_simple(t + dt,   Y + dt*k3,   NPY, NIN);
            Y = Y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        end
    end

    % ---- 2. spike detection  ----
    Vnow = Y([Vidx_PY, Vidx_IN]);
    spiked = (Vprev < 0) & (Vnow >= 0);
    Vprev = Vnow;
    spikedPY = spiked(1:NPY);
    spikedIN = spiked(NPY+1:end);

    if any(spiked)
        ids = find(spiked);
        rasterT = [rasterT; repmat(t, numel(ids),1)];
        rasterID = [rasterID; ids];
    end

    % ---- 3. recurrent AMPA/GABA jumps, full magnitude ----

    % recurrent AMPA (existing) + NMDA (new)
    Y(gAMPA_idx_PY)  = Y(gAMPA_idx_PY)  + 7.5e-3 * (C_PYPY * spikedPY);
    Y(gNMDAs_idx_PY) = Y(gNMDAs_idx_PY) + gNMDA_scale*0.4*7.5e-3 * (C_PYPY * spikedPY); 
    Y(gNMDAf_idx_PY) = Y(gNMDAf_idx_PY) + gNMDA_scale*0.4*7.5e-3 * (C_PYPY * spikedPY); 
    Y(gAMPA_idx_IN)  = Y(gAMPA_idx_IN)  + 2e-3 * (C_PYIN * spikedPY);
    Y(gNMDAs_idx_IN) = Y(gNMDAs_idx_IN) + gNMDA_scale*0.1*2e-3 * (C_PYIN * spikedPY);     
    Y(gNMDAf_idx_IN) = Y(gNMDAf_idx_IN) + gNMDA_scale*0.1*2e-3 * (C_PYIN * spikedPY);     
    % external AMPA (existing) + NMDA (new)
    Y(gAMPA_idx_PY)  = Y(gAMPA_idx_PY)  + g_EXT_PYPY * extSpikePY;
    Y(gNMDAs_idx_PY) = Y(gNMDAs_idx_PY) + gNMDA_scale*0.4*g_EXT_PYPY * extSpikePY;      
    Y(gNMDAf_idx_PY) = Y(gNMDAf_idx_PY) + gNMDA_scale*0.4*g_EXT_PYPY * extSpikePY;       
    Y(gAMPA_idx_IN)  = Y(gAMPA_idx_IN)  + g_EXT_INPY * extSpikeIN;
    Y(gNMDAs_idx_IN) = Y(gNMDAs_idx_IN) + gNMDA_scale*0.1*g_EXT_INPY * extSpikeIN;      
    Y(gNMDAf_idx_IN) = Y(gNMDAf_idx_IN) + gNMDA_scale*0.1*g_EXT_INPY * extSpikeIN;  

    % GABA terminal update 
    if gaba_simulation == "ode"
        [T_X,T_Y,T_c,T_b] = gaba_terminal_step(T_X,T_Y,T_c,T_b,spikedIN,dt, U_gaba,Util_gaba,tau_D,tau_R,lambda_t,Ka,alpha_ca,KP,beta_ca,kplus,kminus,bt,IP);
    elseif gaba_simulation == "simple"
        Y(gGABA_idx_PY) = Y(gGABA_idx_PY) + gGABA_scale*0.8e-3 * (C_INPY * spikedIN);
        Y(gGABA_idx_IN) = Y(gGABA_idx_IN) + gGABA_scale*5e-4   * (C_ININ * spikedIN);
        Y(gGABA_idx_PY) = Y(gGABA_idx_PY) + g_EXT_PYIN * extGabaPY;
        Y(gGABA_idx_IN) = Y(gGABA_idx_IN) + g_EXT_ININ * extGabaIN;
    end

    LFP(k) = mean(Vnow);
    zTrace(k) = mean(Y(3:7:7*NPY));
    gGABATrace(k) = mean(Y(gGABA_idx_PY));
    waitbar(k/nSteps, h, sprintf('Processing %d/%d', k, nSteps))
end

figure; plot(rasterT, rasterID, '.k', 'MarkerSize', 4); xlabel('t [ms]'); ylabel('neuron id');

%% ============================================================
%  HELPER PLOTS
%  ============================================================

Fs = 1000/dt;   % dt in ms -> Fs in Hz. With dt=0.05ms this gives Fs=20000 Hz,
                 % which happens to match the paper's stated 20 kHz LFP sampling rate exactly.

figure;
[s,f,t_spec] = spectrogram(LFP - mean(LFP), 256, 240, 512, Fs);
imagesc(t_spec*1000, f, abs(s).^2);
axis xy;
ylim([0 100]);
xlabel('Time [ms]'); ylabel('Frequency [Hz]');
title('LFP spectrogram'); colorbar;

figure; plot((1:nSteps)*dt, LFP);
xlabel('t [ms]'); ylabel('LFP [mV]');
fprintf('LFP range: [%.4f, %.4f], std: %.6f\n', min(LFP), max(LFP), std(LFP));

binSize = 5;  % ms
edges = 0:binSize:T_total;

figure;
subplot(2,1,1);
histogram(rasterT(rasterID<=NPY), edges);
title('PY population rate'); xlabel('t [ms]'); ylabel('spikes/bin');

subplot(2,1,2);
histogram(rasterT(rasterID>NPY), edges);
title('IN population rate'); xlabel('t [ms]'); ylabel('spikes/bin');

figure;
subplot(2,1,1);
plot((1:nSteps)*dt, zTrace);
xlabel('t [ms]'); ylabel('mean z (PY population)');
title('PY adaptation variable over time');

subplot(2,1,2);
plot((1:nSteps)*dt, gGABATrace);
xlabel('t [ms]'); ylabel('mean gGABA (PY population)');
title('PY inhibitory conductance over time');

function check_rhythm_v2(rasterT, rasterID, NPY, t_start, t_end)
    binSize = 1; % finer resolution -- Nyquist now 500Hz, plenty of headroom
    edges = t_start:binSize:t_end;
    psth = histcounts(rasterT(rasterID<=NPY & rasterT>=t_start & rasterT<=t_end), edges);
    [acor, lags] = xcorr(psth - mean(psth), 'coeff');

    figure;
    plot(lags*binSize, acor); xlim([-100 100]); grid on;
    xlabel('lag [ms]'); ylabel('autocorrelation');
    title(sprintf('[%d-%dms]', t_start, t_end));

    % require peaks to be at least 8ms from lag=0 (above gamma-period range)
    % and at least 8ms apart from each other, to skip trivial near-zero-lag bumps
    [pks, locs] = findpeaks(acor(lags>8), 'MinPeakDistance', 8);
    if ~isempty(locs)
        period = lags(lags>8); period = period(locs(1));
        fprintf('[%d-%dms]: first real peak at lag=%.1fms -> %.1fHz (height=%.2f)\n', ...
            t_start, t_end, period, 1000/period, pks(1));
    end
end

function dY_out = computeDY(t, Y, NPY, NIN, C_INPY, C_ININ, g_GABA_PYIN, g_GABA_ININ, Yrel)
    dY_out = zeros(size(Y));
    for i = 1:NPY
        idx = (i-1)*7 + (1:7);
        dY_out(idx) = PY_neuron_ode(t, Y(idx), [0,0,0,0]);
    end
    off = NPY*7;
    for i = 1:NIN
        idx = off + (i-1)*7 + (1:7);
        dY_out(idx) = IN_neuron_ode(t, Y(idx), [0,0,0,0]);
    end
    gGABA_idx_PY = 7:7:7*NPY;
    gGABA_idx_IN = off+7:7:off+7*NIN;
    dY_out(gGABA_idx_PY) = dY_out(gGABA_idx_PY) + g_GABA_PYIN*(C_INPY*Yrel);
    dY_out(gGABA_idx_IN) = dY_out(gGABA_idx_IN) + g_GABA_ININ*(C_ININ*Yrel);
end

check_rhythm_v2(rasterT, rasterID, NPY, 100, 250);
check_rhythm_v2(rasterT, rasterID, NPY, 350, 500);

function dY = computeDY_simple(t, Y, NPY, NIN)
    dY = zeros(size(Y));
    for i = 1:NPY
        idx = (i-1)*7 + (1:7);
        dY(idx) = PY_neuron_ode(t, Y(idx), [0,0,0,0]);
    end
    off = NPY*7;
    for i = 1:NIN
        idx = off + (i-1)*7 + (1:7);
        dY(idx) = IN_neuron_ode(t, Y(idx), [0,0,0,0]);
    end
end
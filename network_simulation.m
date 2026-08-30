clear;
clc;
close all;

%% ============================================================
%  NETWORK PARAMETERS
%  ============================================================

Nx = 3;
Ny = 3;

N = Nx * Ny;
gGABA_scale = 1;
gNMDA_scale = 1;
integration_method = "R4K"; % R4K, Euler
gaba_simulation = "ode"; % simple, ode

rng(1);               % Fixed seed for random number generation

% External drives
vSTIM_PY = 250;
vSTIM_IN = 0;

% Neuron types
% 1 = PY
% 2 = IN

type = ones(N,1);
type(1:5:N) = 2;

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

p_PY_PY = 0; %0.4;
p_PY_IN = 1; %0.3;
p_IN_PY = 1; %0.6;
p_IN_IN = 0; %0.7;

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

%C = zeros(N,N);

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
beta      = 5;
KP        = 0.4;
gama      = 0.08;

kminus    = 0.95;
kplus     = kminus/0.051;
bt        = 100;
c_rest    = 0.05;
IP        = 0.1102; % beta * c_rest^2 / (c_rest^2 + KP^2);

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

% External drives probability calculation
p_ext_PY = vSTIM_PY * dt / 1000;                 
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




%   fprintf('Kicked %d neurons. Vm values: ', n_kick); disp(Y(Vidx_PY(kick_idx))');
spikedPY = false(NPY,1);
spikedIN = false(NIN,1);

% track state variables
Vm_trace = zeros(nSteps, NPY+NIN); 
GABA_trace = zeros(nSteps, NPY+NIN); 

for k = 1:nSteps
    t = k*dt;
    extSpikePY = rand(NPY,1) < p_ext_PY;
    extSpikeIN = rand(NIN,1) < p_ext_IN;
    extGabaPY = rand(NPY,1) < p_ext_PY;   % same rate as extSpikePY, independent draw
    extGabaIN = rand(NIN,1) < p_ext_IN;   % same rate as extSpikeIN, independent draw

    arrivedAMPA_PY = (C_PYPY * spikedPY) > 0;
    %arrivedGABA_PY = (C_INPY * spikedIN) > 0;
    arrivedAMPA_IN = (C_PYIN * spikedPY) > 0;
    %arrivedGABA_IN = (C_ININ * spikedIN) > 0;

    %arrivedAMPA_PY = C_PYPY * spikedPY
    %arrivedAMPA_IN = C_PYIN * spikedPY

    C_PYPY;
    C_PYIN;
    C_INPY;
    C_ININ;

    if any(spikedPY)
        spikedPY;
    end

    % ---- 1. decay-only integration ----
    if integration_method == "Euler"
        dY = networkODE(t, Y, false(NPY+NIN,1), false(NPY+NIN,1), I_amp_PY, I_amp_IN, NPY, NIN);
        Y = Y + dt*dY;
    elseif integration_method == "R4K"
        if gaba_simulation == "ode"
            Yrel = T_Y;   % freeze this step's release fraction for the RK4 substeps
            gGABA_source_PY = gGABA_scale*0.8e-3 * (C_INPY * Yrel);   % NPYx1, continuous rate
            gGABA_source_IN = gGABA_scale*5e-4   * (C_ININ * Yrel);  
            k1 = computeDY_spikes(t,        Y,           NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt);
            k2 = computeDY_spikes(t + dt/2, Y + dt/2*k1, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt);
            k3 = computeDY_spikes(t + dt/2, Y + dt/2*k2, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt);
            k4 = computeDY_spikes(t + dt,   Y + dt*k3,   NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt);
            Y = Y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        elseif gaba_simulation == "simple"
            k1 = computeDY_spikes(t,        Y,           NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt);
            k2 = computeDY_spikes(t + dt/2, Y + dt/2*k1, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt);
            k3 = computeDY_spikes(t + dt/2, Y + dt/2*k2, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt);
            k4 = computeDY_spikes(t + dt,   Y + dt*k3,   NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt);
            Y = Y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        end
    end

    % ---- 2. spike detection  ----
    Vnow = Y([Vidx_PY, Vidx_IN]);
    
    if (k==1)
        % Force a subset of neurons to spike at t=0, to seed activity
        %kick_fraction = 1;   % kick 10% of PY neurons, e.g.
        %n_kick = round(kick_fraction * NPY);
        %kick_idx = randperm(NPY, n_kick);
        %Y(Vidx_PY(kick_idx)) = 70;
        %Vnow = Y([Vidx_PY, Vidx_IN]);
        %fprintf('DEBUG kick: n_kick=%d, Vnow(kick_idx(1))=%.2f, Vprev(kick_idx(1))=%.2f\n', ...
        %n_kick, Vnow(kick_idx(1)), Vprev(kick_idx(1)));
    end

    if (k==1)
        %kick_fraction_IN = 1;   % set independently — doesn't have to match PY
        %n_kick_IN = round(kick_fraction_IN * NIN);
        %kick_idx_IN = randperm(NIN, n_kick_IN);
        %Y(Vidx_IN(kick_idx_IN)) = 70;
    
        Vnow = Y([Vidx_PY, Vidx_IN]);
    end
    Vm_trace(k,:) = Vnow;

    spiked = (Vprev < 0) & (Vnow >= 0);
    Vprev = Vnow;
    spikedPY = spiked(1:NPY);
    spikedIN = spiked(NPY+1:end);

    if any(spiked)
        ids = find(spiked);
        rasterT = [rasterT; repmat(t, numel(ids),1)];
        rasterID = [rasterID; ids];
    end

    if any(spiked)
        fprintf('t=%.2f: %d neurons spiked\n', t, sum(spiked));
    end

  

    % GABA terminal update 
    if gaba_simulation == "ode"
        [T_X,T_Y,T_c,T_b] = gaba_terminal_step(T_X,T_Y,T_c,T_b,spikedIN,dt, U_gaba,Util_gaba,tau_D,tau_R,lambda_t,Ka,beta,KP,gama,kplus,kminus,bt,IP);
    elseif gaba_simulation == "simple"
        %Y(gGABA_idx_PY) = Y(gGABA_idx_PY) + gGABA_scale*0.8e-3 * (C_INPY * spikedIN);
        %Y(gGABA_idx_IN) = Y(gGABA_idx_IN) + gGABA_scale*5e-4   * (C_ININ * spikedIN);
        %Y(gGABA_idx_PY) = Y(gGABA_idx_PY) + g_EXT_PYIN * extGabaPY;
        %Y(gGABA_idx_IN) = Y(gGABA_idx_IN) + g_EXT_ININ * extGabaIN;
    end

    LFP(k) = mean(Vnow);
    zTrace(k) = mean(Y(3:7:7*NPY));
    gGABATrace(k) = mean(Y(gGABA_idx_PY));
    GABA_trace(k,:) = Y(7:7:end);
    waitbar(k/nSteps, h, sprintf('Processing %d/%d', k, nSteps))

    if k ==1
        % Right after the very first loop iteration (or check at k=2, after the jump from step 1 has applied):
        %fprintf('Max gAMPA among PY neurons after step 1: %.6f\n', max(Y(gAMPA_idx_PY)));
        %fprintf('Number of PY neurons with nonzero gAMPA: %d\n', sum(Y(gAMPA_idx_PY) > 0));
        %fprintf('DEBUG spike-check k=1: sum(spiked)=%d, sum(spikedPY)=%d\n', sum(spiked), sum(spikedPY));
    end

    if k <= 200   % first 10ms at dt=0.05 — plenty to see one spike propagate
        fprintf('k=%3d t=%.2f | spiked=[%s] | arrivedAMPA=[%s] | gAMPA=[%s]\n | Vm=[%s]\n', ...
        k, t, ...
        sprintf('%d',spikedPY), ...
        sprintf('%d',arrivedAMPA_PY), ...
        sprintf('%.4f ',Y(gAMPA_idx_PY)), ...
        sprintf('%.5f',Y(1:7:end)));
    end
    %fprintf('k=%d t=%.2f: V(kick neuron)=%.2f\n', k, t, Y(Vidx_PY(kick_idx(1))));
end

figure; plot(rasterT, rasterID, '.k'); xlabel('t [ms]'); ylabel('neuron id');

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

function dY = computeDY_simple(t, Y, NPY, NIN, spikedPY_flag, spikedIN_flag)
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

function dY = computeDY_spikes(t, Y, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt)
    dY = zeros(size(Y));
    for i = 1:NPY
        idx = (i-1)*7 + (1:7);
        p_i = [0, 0, arrivedAMPA_PY(i), gGABA_source_PY(i), extSpikePY(i), extGabaPY(i), dt];
        dY(idx) = PY_neuron_ode(t, Y(idx), p_i);
    end
    off = NPY*7;
    for i = 1:NIN
        idx = off + (i-1)*7 + (1:7);
        p_i = [0, 0, arrivedAMPA_IN(i), gGABA_source_IN(i), extSpikeIN(i), extGabaIN(i), dt];
        dY(idx) = IN_neuron_ode(t, Y(idx), p_i);
    end
end

figure;
offset_step = 100;  % mV spacing between neurons, adjust to taste
hold on;
timeAxis = (1:nSteps)*dt;
for i = 1:(NPY+NIN)
    plot(timeAxis, Vm_trace(:,i) + (i-1)*offset_step, 'k');
end
% mark PY/IN boundary
yline((NPY-0.5)*offset_step, '--r');
xlabel('t [ms]');
ylabel('Neuron id (offset for visibility)');
title('Vm per neuron over time (PY below dashed line, IN above)');
yticks((0:(NPY+NIN-1))*offset_step);
yticklabels(1:(NPY+NIN));
hold off;

figure;
offset_step = 100;  % mV spacing between neurons, adjust to taste
hold on;
timeAxis = (1:nSteps)*dt;
for i = 1:(NPY+NIN)
    plot(timeAxis, GABA_trace(:,i) + (i-1)*offset_step, 'k');
end
% mark PY/IN boundary
yline((NPY-0.5)*offset_step, '--r');
xlabel('t [ms]');
ylabel('Neuron id (offset for visibility)');
title('gGABA per neuron over time (PY below dashed line, IN above)');
yticks((0:(NPY+NIN-1))*offset_step);
yticklabels(1:(NPY+NIN));
hold off;
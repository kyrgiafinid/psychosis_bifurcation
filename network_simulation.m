clear;
clc;
close all;

%% ============================================================
%  NETWORK PARAMETERS
%  ============================================================

Nx = 30;
Ny = 30;

N = Nx * Ny;
gGABA_scale = 1;
gNMDA_scale = 1;
integration_method = "R4K"; % R4K, Euler
gaba_simulation = "ode"; % simple, ode

rng(1);               % Fixed seed for random number generation

% External drives
vSTIM_PY = 250;
vSTIM_IN = 500;

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
        dx = min(dx, Nx - dx);           % wraparound

        dy = abs(y(pre) - y(post));
        dy = min(dy, Ny - dy);           % wraparound

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

nStatesIN = 11;
nStatesPY = 7;

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
      0;            % gGABA   - no synapse active yet
      1;            % X       - fraction of GABA resource "ready to release"
      0;            % Y       - fraction of GABA resource "actively being released" right now
      c_rest;       % c       - residual presynaptic calcium concentration
      0];           % b       - calcium bound to parvalbumin


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
T_total = 1000;
nSteps = round(T_total/dt);

Y = [repmat(Y_PY, NPY, 1); repmat(Y_IN, NIN, 1)];
Vidx_PY = 1:nStatesPY:nStatesPY*NPY;
Vidx_IN = nStatesPY*NPY+1 : nStatesIN : nStatesPY*NPY+nStatesIN*NIN;

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


Xidx_IN = 7*NPY+8 : nStatesIN : 7*NPY+nStatesIN*NIN;
Yidx_IN = 7*NPY+9 : nStatesIN : 7*NPY+nStatesIN*NIN;
cidx_IN = 7*NPY+10: nStatesIN : 7*NPY+nStatesIN*NIN;
bidx_IN = 7*NPY+11: nStatesIN : 7*NPY+nStatesIN*NIN;
gAMPA_idx_IN = 7*NPY+4:7:7*NPY+7*NIN;
gGABA_idx_PY = 7:7:7*NPY;
gGABA_idx_IN = 7*NPY+7 : nStatesIN : 7*NPY+nStatesIN*NIN;
gNMDAs_idx_PY = 5:7:7*NPY;  gNMDAf_idx_PY = 6:7:7*NPY;
gNMDAs_idx_IN = 7*NPY+5:7:7*NPY+7*NIN;  gNMDAf_idx_IN = 7*NPY+6:7:7*NPY+7*NIN;
bt_vec = bt * ones(NIN,1); 

% testing plots 
zTrace = zeros(nSteps,1);
gGABATrace = zeros(nSteps,1);
gGABAsourceTrace = zeros(nSteps,1);
YrelMeanTrace = zeros(nSteps,1);

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
    arrivedAMPA_IN = (C_PYIN * spikedPY) > 0;

    arrivedAMPA_PY = C_PYPY * spikedPY;
    arrivedAMPA_IN = C_PYIN * spikedPY;
    
    % -------GABA release terminal constants (Tsodyks-Markram + calcium/PV)
    lambda_t = 0.03e-3;   % async-release rate constant
    Ka = 0.2;              % calcium affinity of async release
    
    c_prev = Y(cidx_IN);   % last step's calcium level, per IN neuron
    lam = lambda_t .* (c_prev.^4) ./ (c_prev.^4 + Ka.^4);
    asyncEvent = rand(NIN,1) < (1 - exp(-lam.*dt));
    bt_vec = bt * ones(NIN,1);

    Yrel = Y(Yidx_IN);  
    gGABA_source_PY = gGABA_scale*0.8e-3 * (C_INPY * Yrel);   % NPYx1, continuous rate
    gGABA_source_IN = gGABA_scale*5e-4   * (C_ININ * Yrel);  
    k1 = computeDY_spikes_gaba(t,        Y,           NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, gNMDA_scale,bt_vec, spikedIN,asyncEvent);
    k2 = computeDY_spikes_gaba(t + dt/2, Y + dt/2*k1, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, gNMDA_scale,bt_vec, spikedIN,asyncEvent);
    k3 = computeDY_spikes_gaba(t + dt/2, Y + dt/2*k2, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, gNMDA_scale,bt_vec, spikedIN,asyncEvent);
    k4 = computeDY_spikes_gaba(t + dt,   Y + dt*k3,   NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, gNMDA_scale,bt_vec, spikedIN,asyncEvent);
    Y = Y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    % ---- 2. spike detection  ----
    Vnow = Y([Vidx_PY, Vidx_IN]);
    
    Vm_trace(k,:) = Vnow;
    gGABAsourceTrace(k) = mean(gGABA_source_PY);
    YrelMeanTrace(k) = mean(Yrel);

    spiked = (Vprev < 0) & (Vnow >= 0);
    Vprev = Vnow;
    spikedPY = spiked(1:NPY);
    spikedIN = spiked(NPY+1:end);

    if any(spiked)
        ids = find(spiked);
        rasterT = [rasterT; repmat(t, numel(ids),1)];
        rasterID = [rasterID; ids];
    end

    gGABA_idx_PY = 7:7:7*NPY;                      
    gGABA_idx_IN = 7*NPY+7 : 11 : 7*NPY+11*NIN;  
    LFP(k) = mean(Vnow);
    zTrace(k) = mean(Y(3:7:7*NPY));
    gGABATrace(k) = mean(Y(gGABA_idx_PY));
    GABA_trace(k,:) = [Y(gGABA_idx_PY)', Y(gGABA_idx_IN)'];
    waitbar(k/nSteps, h, sprintf('Processing %d/%d', k, nSteps))

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

figure;
subplot(2,1,1); plot((1:nSteps)*dt, gGABAsourceTrace); ylabel('mean gGABA\_source\_PY'); title('Source term over time');
subplot(2,1,2); plot((1:nSteps)*dt, YrelMeanTrace); ylabel('mean Yrel (IN release var)'); title('IN Y state over time');
title('GABA and Yrel');

function check_rhythm_v2(rasterT, rasterID, NPY, t_start, t_end)
    binSize = 1; % finer resolution -- Nyquist now 500Hz, plenty of headroom
    edges = t_start:binSize:t_end;
    psth = histcounts(rasterT(rasterID<=NPY & rasterT>=t_start & rasterT<=t_end), edges);
    [acor, lags] = xcorr(psth - mean(psth), 'coeff');

    figure;
    plot(lags*binSize, acor); xlim([-100 100]); grid on;
    xlabel('lag [ms]'); ylabel('autocorrelation');
    title(sprintf('[%d-%dms]', t_start, t_end));

    % require peaks to be at least 8ms from lag=0 (above gama-period range)
    % and at least 8ms apart from each other, to skip trivial near-zero-lag bumps
    [pks, locs] = findpeaks(acor(lags>8), 'MinPeakDistance', 8);
    if ~isempty(locs)
        period = lags(lags>8); period = period(locs(1));
        fprintf('[%d-%dms]: first real peak at lag=%.1fms -> %.1fHz (height=%.2f)\n', ...
            t_start, t_end, period, 1000/period, pks(1));
    end
end



check_rhythm_v2(rasterT, rasterID, NPY, 100, 250);
check_rhythm_v2(rasterT, rasterID, NPY, 350, 500);

if (N <= 25)
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
end

if (N <= 25)
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
end

function dY = computeDY_spikes_gaba(t, Y, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, gNMDA_scale, bt_vec, spikedIN, asyncEvent)
    dY = zeros(size(Y));

    % PY: 7 states each, 11-element p (last 3 unused, kept for layout parity)
    for i = 1:NPY
        idx = (i-1)*7 + (1:7);
        p_i = [0, 0, arrivedAMPA_PY(i), gGABA_source_PY(i), extSpikePY(i), extGabaPY(i), dt, gNMDA_scale, 0, 0, 0];
        dY(idx) = PY_neuron_ode(t, Y(idx), p_i);
    end

    % IN: 11 states each, full 11-element p including bt/spike_own_IN/asyncEvent_own
    off = NPY*7;
    for i = 1:NIN
        idx = off + (i-1)*11 + (1:11);
        p_i = [0, 0, arrivedAMPA_IN(i), gGABA_source_IN(i), extSpikeIN(i), extGabaIN(i), dt, gNMDA_scale, bt_vec(i), spikedIN(i), asyncEvent(i)];
        dY(idx) = IN_neuron_ode(t, Y(idx), p_i);
    end
end

fprintf('t=%.2f: mean gGABA_source_PY = %.6f\n', t, mean(gGABA_source_PY));

figure;
subplot(2,1,1); plot((1:nSteps)*dt, gGABAsourceTrace); ylabel('mean gGABA\_source\_PY'); title('Source term over time');
subplot(2,1,2); plot((1:nSteps)*dt, YrelMeanTrace); ylabel('mean Yrel (IN release var)'); title('IN Y state over time');
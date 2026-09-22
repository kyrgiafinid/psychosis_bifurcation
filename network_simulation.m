clear;
clc;
close all;

%% ============================================================
%  NETWORK PARAMETERS
%  ============================================================

Nx = 30;
Ny = 30;

N = Nx * Ny;

rng(1);               % Fixed seed for random number generation

% External drives
vSTIM_PY = 200;
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

figure('Name','Lattice');

scatter(x(PY),y(PY),50,'b','filled');
hold on;

scatter(x(IN),y(IN),70,'r','filled');

xlabel('x');
ylabel('y');

title('Neuronal Lattice');

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

figure('Name','Connectivity matrix');

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

g_GABA_IN = 0.5;
g_GABA_PY = 0.8;

% Parameters for Tsodyks-Markram + calcium/PV GABA release parameters (paper eqs. 13-18)
lambda_t = 0.03;
Ka = 0.2;
bt = 100;
c_rest = 0.05;
bt_vec = bt * ones(NIN,1);

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

%% Initializations
% ---- reorder connectivity to match Y's [PY block; IN block] layout ----
C_PYPY = C(PY,PY);   % PY - PY connections
C_INPY = C(PY,IN);   % PY - IN connections
C_PYIN = C(IN,PY);   % IN - PY connections
C_ININ = C(IN,IN);   % IN - IN connections

N_in_PYPY = max(sum(C_PYPY, 2), 1);   % PY inputs per PY
N_in_PYIN = max(sum(C_PYIN, 2), 1);   % PY inputs per IN
N_in_INPY = max(sum(C_INPY, 2), 1);   % IN inputs per PY
N_in_ININ = max(sum(C_ININ, 2), 1);   % IN inputs per IN

I_amp_PY = zeros(NPY,1);
I_amp_IN = zeros(NIN,1);

dt = 0.05;          
T_total = 1000;
nSteps = round(T_total/dt);

%% preallocations
Y = [repmat(Y_PY, NPY, 1); repmat(Y_IN, NIN, 1)];
Vidx_PY = 1:nStatesPY:nStatesPY*NPY;
Vidx_IN = nStatesPY*NPY+1 : nStatesIN : nStatesPY*NPY+nStatesIN*NIN;

Vprev = Y([Vidx_PY, Vidx_IN]);   % for threshold-crossing detection

spikeAMPA = false(NPY+NIN,1);
spikeGABA = false(NPY+NIN,1);

maxSpikesExpected = nSteps * 50; 
rasterT = zeros(maxSpikesExpected, 1);
rasterID = zeros(maxSpikesExpected, 1);
rasterCount = 0;

spikedPY = false(NPY,1);
spikedIN = false(NIN,1);

totalAsyncEventCount = 0;

% External drives probability calculation
p_ext_PY = 1 - exp(-vSTIM_PY * dt / 1000);
p_ext_IN = 1 - exp(-vSTIM_IN * dt / 1000);

g_EXT_PYPY = 0.25;    % paper's value, PY's external excitatory conductance
g_EXT_INPY = 0.003;   % paper's value, IN's external excitatory conductance
g_EXT_PYIN = 0.025;   % external inhibitory conductance onto PY (g~EXT^{PY<-IN})
g_EXT_ININ = 1e-4;    % external inhibitory conductance onto IN (g~EXT^{IN<-IN}


% trace states
Xidx_IN = 7*NPY+8 : nStatesIN : 7*NPY+nStatesIN*NIN;
Yidx_IN = 7*NPY+9 : nStatesIN : 7*NPY+nStatesIN*NIN;
cidx_IN = 7*NPY+10: nStatesIN : 7*NPY+nStatesIN*NIN;
bidx_IN = 7*NPY+11: nStatesIN : 7*NPY+nStatesIN*NIN;
gGABA_idx_PY = 7:7:7*NPY;
gGABA_idx_IN = 7*NPY+7 : nStatesIN : 7*NPY+nStatesIN*NIN;
gNMDAs_idx_PY = 5:7:7*NPY;  gNMDAf_idx_PY = 6:7:7*NPY;

%% preparation for state and parameter tracking
LFP = zeros(nSteps,1);

gAMPA_idx_PY  = 4:7:7*NPY;
gAMPA_idx_IN  = 7*NPY+4 : nStatesIN : 7*NPY+nStatesIN*NIN;
gGABA_all_trace = zeros(nSteps,1);
Vm_trace = zeros(nSteps, NPY+NIN);

%h = waitbar(0, 'Processing...');
PY_idx_all = arrayfun(@(i) (i-1)*7+(1:7), 1:NPY, 'UniformOutput', false);
IN_idx_all = arrayfun(@(i) NPY*7+(i-1)*11+(1:11), 1:NIN, 'UniformOutput', false);

%vectorization arrays
idxPY.V=1:7:7*NPY; idxPY.w=2:7:7*NPY; idxPY.z=3:7:7*NPY;
idxPY.gAMPA=4:7:7*NPY; idxPY.gNMDAs=5:7:7*NPY; idxPY.gNMDAf=6:7:7*NPY; idxPY.gGABA=7:7:7*NPY;
off=7*NPY;
idxIN.V=off+1:11:off+11*NIN; idxIN.h=off+2:11:off+11*NIN; idxIN.n=off+3:11:off+11*NIN;
idxIN.gAMPA=off+4:11:off+11*NIN; idxIN.gNMDAs=off+5:11:off+11*NIN; idxIN.gNMDAf=off+6:11:off+11*NIN;
idxIN.gGABA=off+7:11:off+11*NIN; idxIN.X=off+8:11:off+11*NIN; idxIN.Y=off+9:11:off+11*NIN;
idxIN.c=off+10:11:off+11*NIN; idxIN.b=off+11:11:off+11*NIN;
for k = 1:nSteps
    
    % -------------1. Spiking probabilities -------------
    t = k*dt;
    extSpikePY = rand(NPY,1) < p_ext_PY;
    extSpikeIN = rand(NIN,1) < p_ext_IN;
    extGabaPY = rand(NPY,1) < p_ext_PY;   
    extGabaIN = rand(NIN,1) < p_ext_IN;   
        
    Yrel = Y(Yidx_IN);  
    arrivedAMPA_PY = (C_PYPY * spikedPY);
    arrivedAMPA_IN = (C_PYIN * spikedPY);
    gGABA_source_PY = g_GABA_PY * (C_INPY * Yrel);
    gGABA_source_IN = g_GABA_IN * (C_ININ * Yrel);
        
    % probability of asynchronus gaba spiking   
    c_prev = Y(cidx_IN);  
    lam = lambda_t .* (c_prev.^4) ./ (c_prev.^4 + Ka.^4);
    asyncEvent = 0.01 * rand(NIN,1) < (1 - exp(-lam.*dt));
    
    % -------------2. R4K implementation ----------------
    k1 = compute_step_vec(t,        Y,           NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, bt_vec, spikedIN, asyncEvent, idxPY, idxIN);
    k2 = compute_step_vec(t + dt/2, Y + dt/2*k1, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, bt_vec, spikedIN, asyncEvent, idxPY, idxIN);
    k3 = compute_step_vec(t + dt/2, Y + dt/2*k2, NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, bt_vec, spikedIN, asyncEvent, idxPY, idxIN);
    k4 = compute_step_vec(t + dt,   Y + dt*k3,   NPY, NIN, arrivedAMPA_PY, arrivedAMPA_IN, gGABA_source_PY, gGABA_source_IN, extSpikePY, extGabaPY, extSpikeIN, extGabaIN, dt, bt_vec, spikedIN, asyncEvent, idxPY, idxIN);
    Y = Y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    % -------------3. Spike detection  -------------------
    Vnow = Y([Vidx_PY, Vidx_IN]); 
    Vm_trace(k,:) = Vnow;

    spiked = (Vprev < 0) & (Vnow >= 0);
    Vprev = Vnow;

    spikedPY = spiked(1:NPY);
    spikedIN = spiked(NPY+1:end);

    if any(spiked)
        ids = find(spiked);
        n = numel(ids);
        rasterT(rasterCount+1 : rasterCount+n) = t;
        rasterID(rasterCount+1 : rasterCount+n) = ids;
        rasterCount = rasterCount + n;
    end
    
    % -------------3. Track variables  -------------------
    gGABA_all_trace(k) = mean([Y(gGABA_idx_PY); Y(gGABA_idx_IN)]);
    LFP(k) = mean(Vnow);

    %waitbar(k/nSteps, h, sprintf('Processing %d/%d', k, nSteps))
    if mod(k, 200) == 0
        fprintf('Processing %d/%d\n', k, nSteps);
    end
end

figure; plot(rasterT, rasterID, '.k'); xlabel('t [ms]'); ylabel('neuron id');

%% ============================================================
%  HELPER PLOTS
%  ============================================================

Fs = 1000/dt;  
rasterT = rasterT(1:rasterCount);
rasterID = rasterID(1:rasterCount);

figure('Name', 'mean LPF af all neurons');
plot((1:nSteps)*dt, LFP);
xlabel('t [ms]'); ylabel('LFP [mV]');
fprintf('LFP range: [%.4f, %.4f], std: %.6f\n', min(LFP), max(LFP), std(LFP));

binSize = 5;  % ms
edges = 0:binSize:T_total;

figure('Name', 'histogram of firing neurons');
subplot(2,1,1);
histogram(rasterT(rasterID<=NPY), edges);
title('PY population rate'); xlabel('t [ms]'); ylabel('spikes/bin');
subplot(2,1,2);
histogram(rasterT(rasterID>NPY), edges);
title('IN population rate'); xlabel('t [ms]'); ylabel('spikes/bin');

figure('Name', ' mean gGABA of all neurons in the network');
plot((1:nSteps)*dt, gGABA_all_trace);
xlabel('t [ms]');
ylabel('mean gGABA (all 900 neurons)');
title('Inhibitory conductance over time, PY+IN combined');


if (N <= 25)
    figure('Name', 'plot membrane potential of all neurons');
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


window_len = 2000;   % samples = 2000/20000*1000 = 100ms of data per window -> ~10Hz resolution
skip_ms = 100;   % exclude the initial transient (visible in your plot, roughly 0-100ms)
skip_idx = round(skip_ms/dt);
[s,f,t_spec] = spectrogram(LFP(skip_idx:end) - mean(LFP(skip_idx:end)), window_len, round(window_len*0.9), 2048, Fs);

figure('Name','Spectrogram');
imagesc(t_spec*1000 + skip_ms, f, abs(s).^2);
axis xy;
ylim([0 100]);
xlabel('Time [ms]'); ylabel('Frequency [Hz]');
title('LFP spectrogram (Fig. 1C equivalent)');
colorbar;

avg_power = mean(abs(s).^2, 2);
[peak_val, peak_idx] = max(avg_power);
fprintf('Dominant frequency (max time-averaged power): %.1f Hz (power=%.2e)\n', f(peak_idx), peak_val);

rate_PY = sum(rasterID <= NPY) / NPY / (T_total/1000);
rate_IN = sum(rasterID >  NPY) / NIN / (T_total/1000);

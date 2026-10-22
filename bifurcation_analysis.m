%% Bifurcation analysis of wilson Cowan model
%
% System:
%   dE/dt = (-E + (1 - p.rE*E) / (1 + exp(-p.k_E * (xE - p.x0_E))) - 1/(1 + exp(p.k_E * p.x0_E)))/p.tau_E
%   dI/dt = (-I + (1 - p.rI*I) / (1 + exp(-p.k_I * (xI - p.x0_I))) - 1/(1 + exp(p.k_I * p.x0_I)))/p.tau_I

%% Step 1: Define parameters and initial equilibrium
clear;
close all;
p = struct();
p.k_E = 0.1603;         p.x0_E = 18;   
p.k_I = 2.1669;         p.x0_I = 0.949;
p.rE  = 0.01820;        p.rI   = 0.05880;
p.tau_E = 2;            p.tau_I = 8;
p.P_E = 16;             p.P_I  = 0.261;      
p.w_EE = 85;            p.w_EI = 85;

% initial values for w_EI and w_II
p.w_IE = 12.5;            p.w_II = 7.5;

% Find initial equilibrium & and verify instability by chacking eigenvalues
% of the jacobian
F = @(x) [
    (-x(1) + (1 - p.rE*x(1)) * (1 / (1 + exp(-p.k_E * (p.w_EE*x(1) - p.w_EI*x(2) + p.P_E - p.x0_E))) - 1/(1 + exp(p.k_E * p.x0_E))))/p.tau_E
    (-x(2) + (1 - p.rI*x(2)) * (1 / (1 + exp(-p.k_I * (p.w_IE*x(1) - p.w_II*x(2) + p.P_I - p.x0_I))) - 1/(1 + exp(p.k_I * p.x0_I))))/p.tau_I
];

opts_fs = optimoptions('fsolve', 'Display', 'off');
x0 = fsolve(F, [0.15; 0.05], opts_fs);

fprintf('======== 1.Find Equillibria ========\n\n');
fprintf('  E = %.6f\n  I = %.6f\n\n', x0(1), x0(2));

SE0  = 1/(1 + exp(-p.k_E*(p.w_EE*x0(1) - p.w_EI*x0(2) + p.P_E - p.x0_E)));    % E sigmoid at the equilibrium
SI0  = 1/(1 + exp(-p.k_I*(p.w_IE*x0(1) - p.w_II*x0(2) + p.P_I - p.x0_I)));    % I sigmoid at the equilibrium
cE0  = 1/(1 + exp(p.k_E*p.x0_E));                                             % zero-correction constants
cI0  = 1/(1 + exp(p.k_I*p.x0_I));
dSE0 = p.k_E*SE0*(1 - SE0);                                                   % their slopes
dSI0 = p.k_I*SI0*(1 - SI0);
J0 = [(-1 - p.rE*(SE0-cE0) + (1 - p.rE*x0(1))*dSE0*p.w_EE)/p.tau_E,   -(1 - p.rE*x0(1))*dSE0*p.w_EI/p.tau_E;
        (1 - p.rI*x0(2))*dSI0*p.w_IE/p.tau_I,                    (-1 - p.rI*(SI0-cI0) - (1 - p.rI*x0(2))*dSI0*p.w_II)/p.tau_I];
eigs0 = eig(J0);
fprintf('Eigenvalues at start: %.4f + %.4fi,  %.4f + %.4fi\n\n', ...
    real(eigs0(1)), imag(eigs0(1)), real(eigs0(2)), imag(eigs0(2)));

%% Step 2: Equilibrium continuation
function out = continuation_equations
    out = {[]; @fun_eval; @jacobian; @jacobianp; @hessians; []; []; []; []};
end

function dydt = fun_eval(~,x,w_EE,w_EI,w_IE,w_II,P_E,P_I,k_E,x0_E,k_I,x0_I,rE,rI,tau_E,tau_I)
    dydt = [(-x(1) + (1 - rE*x(1)) * (1 / (1 + exp(-k_E * (w_EE*x(1) - w_EI*x(2) + P_E - x0_E))) - 1/(1 + exp(k_E * x0_E))))/tau_E
            (-x(2) + (1 - rI*x(2)) * (1 / (1 + exp(-k_I * (w_IE*x(1) - w_II*x(2) + P_I - x0_I))) - 1/(1 + exp(k_I * x0_I))))/tau_I];
end

opt = contset;
opt = contset(opt, 'Singularities', 1);
opt = contset(opt, 'MaxNumPoints', 10000);
opt = contset(opt, 'MaxStepsize', 0.2);

ap1 = 4;
par = [p.w_EE; p.w_EI; p.w_IE; p.w_II; p.P_E; p.P_I; p.k_E; p.x0_E; p.k_I; p.x0_I; p.rE; p.rI; p.tau_E; p.tau_I];
opt = contset(opt, 'Backward', 1);
[x1, v1]     = init_EP_EP(@continuation_equations, x0, par, ap1);
[x,v,s,h,f]  = cont(@equilibrium, x1, [], opt);

figure('Name','Equilibrium branch (w_II)', 'NumberTitle','off');
cpl(x,v,s,[3,1]);
xlabel('w_{II}', 'FontSize', 12);
ylabel('E', 'FontSize', 12);
title('Equilibrium branch vs w_{II} (decreasing)', 'FontSize', 13);
xlim([0 10]);
hold on;
markerSpecs = struct('H',{{'o','r'}}, 'LP',{{'s','m'}}, 'BT',{{'^','c'}}, 'GH',{{'d','y'}});
legendH = []; legendLabels = {};
for k = 1:length(s)
    lab = s(k).label;
    if any(strcmp(lab, {'00','99'})), continue; end   % skip plain start/end markers
    wIIval = x(3, s(k).index);
    Eval   = x(1, s(k).index);
    if wIIval < 0 || wIIval > 10, continue; end        % stay within the plotted window
    if isfield(markerSpecs, lab)
        spec = markerSpecs.(lab);
    else
        spec = {'x','k'};
    end
    hpt = plot(wIIval, Eval, spec{1}, 'MarkerFaceColor', spec{2}, 'MarkerEdgeColor','k', 'MarkerSize', 9);
    text(wIIval+0.15, Eval, lab, 'FontWeight','bold');
    if ~any(strcmp(legendLabels, lab))
        legendH(end+1) = hpt; legendLabels{end+1} = lab; %#ok<SAGROW>
    end
end
if ~isempty(legendH)
    legend(legendH, legendLabels, 'Location','best');
end

idx = s(2).index;
for i = max(1,idx-3):min(size(x,2),idx+3)
    par(4) = x(3,i);
    Jtest = jacobian(0, x(1:2,i), par(1),par(2),par(3),par(4),par(5),par(6),par(7),par(8),par(9),par(10),par(11),par(12),par(13),par(14));
    fprintf('i=%d  w_II=%.4f  trace=%+.5f\n', i, x(3,i), trace(Jtest));
end

%% Step 3: Hopf continuation
x1 = x(1:2, s(2).index);
par(ap1) = x(end, s(2).index);
[x0,v0] = init_H_H(@continuation_equations, x1, par, [2,4]);

opt = contset(opt, 'Backward', 0);
[x2f,v2f,s2f,h2f,f2f] = cont(@hopf, x0, v0, opt);

opt = contset(opt, 'Backward', 1);
[x2b,v2b,s2b,h2b,f2b] = cont(@hopf, x0, v0, opt);

k = 2;                               % the index where label is 'H'
x1 = x(1:2, s(k).index);
par(ap1) = x(end, s(k).index);
[x0,v0] = init_H_H(@continuation_equations, x1, par, [2,4]);   % w_EI, w_II

opt = contset(opt, 'Backward', 0);
[x2f,v2f,s2f] = cont(@hopf, x0, v0, opt);
opt = contset(opt, 'Backward', 1);
[x2b,v2b,s2b] = cont(@hopf, x0, v0, opt);

% combine both branches of the curve into one ordered list
wEI_c = [fliplr(x2b(3,:)) x2f(3,:)];
wII_c = [fliplr(x2b(4,:)) x2f(4,:)];
Eq_c  = [fliplr(x2b(1:2,:)') ; x2f(1:2,:)'];

delta = 2;  % small perpendicular offset in w_II
side_in  = nan(1,length(wEI_c));
side_out = nan(1,length(wEI_c));
for i = 1:length(wEI_c)
    seed = Eq_c(i,:)';
    side_in(i)  = side_class(wEI_c(i), wII_c(i)-delta, seed, par, opts_fs);
    side_out(i) = side_class(wEI_c(i), wII_c(i)+delta, seed, par, opts_fs);
end

wEI_c = [fliplr(x2b(3,:)) x2f(3,:)];
wII_c = [fliplr(x2b(4,:)) x2f(4,:)];
Eq_c  = [fliplr(x2b(1:2,:)') ; x2f(1:2,:)'];

below_pts = []; above_pts = [];
for i = 1:length(wEI_c)
    seed = Eq_c(i,:)';
    s_below = side_class(wEI_c(i), wII_c(i)-delta, seed, par, opts_fs);
    s_above = side_class(wEI_c(i), wII_c(i)+delta, seed, par, opts_fs);
    if ~isnan(s_below), below_pts = [below_pts; wEI_c(i), wII_c(i)-delta, s_below]; end
    if ~isnan(s_above), above_pts = [above_pts; wEI_c(i), wII_c(i)+delta, s_above]; end
end

seed = x0(1:2);
sign_below = side_class(x0(3), x0(4)-0.6, seed, par, opts_fs);
sign_above = side_class(x0(3), x0(4)+0.6, seed, par, opts_fs);

if sign_below < 0
    text(x0(3)-8, x0(4)-2, 'stable', 'FontSize', 12, 'FontWeight','bold', 'Color',[0 0.3 0.6]);
    text(x0(3)+3, x0(4)+1.5, 'oscillatory', 'FontSize', 12, 'FontWeight','bold', 'Color',[0.7 0.2 0]);
else
    text(x0(3)-8, x0(4)-2, 'oscillatory', 'FontSize', 12, 'FontWeight','bold', 'Color',[0.7 0.2 0]);
    text(x0(3)+3, x0(4)+1.5, 'stable', 'FontSize', 12, 'FontWeight','bold', 'Color',[0 0.3 0.6]);
end

fprintf('below curve: trace sign = %+.0f (consistent in %d%% of samples)\n', mode(below_pts(:,3)), round(100*mean(below_pts(:,3)==mode(below_pts(:,3)))));
fprintf('above curve: trace sign = %+.0f (consistent in %d%% of samples)\n', mode(above_pts(:,3)), round(100*mean(above_pts(:,3)==mode(above_pts(:,3)))));

ap1 = 4;    % w_II
opt = contset(opt, 'Backward', 1);
opt = contset(opt, 'MaxNumPoints', 20000); 
opt = contset(opt, 'MaxStepsize', 0.1);  
plot_1d_bifurcation(80, par, ap1, opt, 'Equilibrium branch vs w_{II} (w_{EI}=80)');
plot_1d_bifurcation(90, par, ap1, opt, 'Equilibrium branch vs w_{II} (w_{EI}=90)');

below_label = cell(1,length(wEI_c));
above_label = cell(1,length(wEI_c));
for i = 1:length(wEI_c)
    seed = Eq_c(i,:)';
    below_label{i} = classify_stability(wEI_c(i), wII_c(i)-delta, seed, par, opts_fs);
    above_label{i} = classify_stability(wEI_c(i), wII_c(i)+delta, seed, par, opts_fs);
end

figure('Name','Two-parameter bifurcation diagram', 'NumberTitle','off'); hold on;
isStableBelow = strcmp(below_label, 'stable');
isOscAbove    = strcmp(above_label, 'oscillatory');

h1 = plot(wEI_c(isStableBelow), wII_c(isStableBelow)-delta, '.', 'Color',[0.75 0.85 1], 'MarkerSize',10);
plot(wEI_c(~isStableBelow), wII_c(~isStableBelow)-delta, '.', 'Color',[1 0.8 0.75], 'MarkerSize',10);
h2 = plot(wEI_c(isOscAbove), wII_c(isOscAbove)+delta, '.', 'Color',[1 0.8 0.75], 'MarkerSize',10);
plot(wEI_c(~isOscAbove), wII_c(~isOscAbove)+delta, '.', 'Color',[0.75 0.85 1], 'MarkerSize',10);

h3 = plot(wEI_c, wII_c, 'k-', 'LineWidth', 2);
h4 = plot(x0(3), x0(4), 'ko', 'MarkerFaceColor','r', 'MarkerSize',7);

xlabel('w_{EI}'); ylabel('w_{II}');
title('Two-parameter bifurcation diagram (w_{EI}, w_{II})');
legend([h1 h2 h3 h4], {'stable','oscillatory','Hopf curve','start point'}, 'Location','southeast');
grid on; box on;

function label = classify_stability(wEI, wII, seed, par, opts_fs)
    par(2) = wEI; par(4) = wII;
    F = @(x) fun_eval(0, x, par(1),par(2),par(3),par(4),par(5),par(6),par(7),par(8),par(9),par(10),par(11),par(12),par(13),par(14));
    xs = fsolve(F, seed, opts_fs);
    J = jacobian(0, xs, par(1),par(2),par(3),par(4),par(5),par(6),par(7),par(8),par(9),par(10),par(11),par(12),par(13),par(14));
    if trace(J) < 0
        label = 'stable';
    else
        label = 'oscillatory';
    end
end

function plot_1d_bifurcation(wEI_fixed, par, ap1, opt, titleStr)
    par(2) = wEI_fixed;
    par(4) = 7.5;
    F = @(x) fun_eval(0, x, par(1),par(2),par(3),par(4),par(5),par(6),par(7),par(8),par(9),par(10),par(11),par(12),par(13),par(14));
    opts_fs = optimoptions('fsolve','Display','off');
    x0 = fsolve(F, [0.15; 0.05], opts_fs);
    
    [x1,v1] = init_EP_EP(@continuation_equations, x0, par, ap1);
    [x,v,s,h,f] = cont(@equilibrium, x1, [], opt);

    % classify each point as stable/unstable via trace sign
    npts = size(x,2);
    trVals = nan(1,npts);
    for i = 1:npts
        parLoc = par; parLoc(ap1) = x(end,i);
        J = jacobian(0, x(1:2,i), parLoc(1),parLoc(2),parLoc(3),parLoc(4),parLoc(5),parLoc(6),parLoc(7),parLoc(8),parLoc(9),parLoc(10),parLoc(11),parLoc(12),parLoc(13),parLoc(14));
        trVals(i) = trace(J);
    end
    stableMask = trVals < 0;

    % restrict to the meaningful w_II window
    mask = x(end,:) >= 0 & x(end,:) <= 10;

    titleStr = sprintf('Equilibrium branch vs w_{II} (w_{EI}=%g)', wEI_fixed);
    figure('Name', titleStr, 'NumberTitle','off'); hold on;
    plot(x(end, mask & stableMask),  x(1, mask & stableMask),  'b.', 'MarkerSize', 8);
    plot(x(end, mask & ~stableMask), x(1, mask & ~stableMask), 'r.', 'MarkerSize', 8);

    for k = 1:length(s)
        if strcmp(s(k).label, 'H') && x(end,s(k).index) >= 0 && x(end,s(k).index) <= 10
            plot(x(end, s(k).index), x(1, s(k).index), 'ko', 'MarkerFaceColor','g', 'MarkerSize', 8);
        end
        if strcmp(s(k).label, 'LP') && x(end,s(k).index) >= 0 && x(end,s(k).index) <= 10
            plot(x(end, s(k).index), x(1, s(k).index), 'ks', 'MarkerFaceColor','m', 'MarkerSize', 8);
        end
    end

    xlim([0 8]);
    xlabel('w_{II}');
    ylabel('E (equilibrium)');
    title(titleStr);
    legend({'stable','unstable','Hopf','LP'}, 'Location','best');
    grid on; box on;
end

function out = side_class(wEI, wII, seed, par, opts_fs)
    par(2) = wEI; par(4) = wII;
    F = @(x) fun_eval(0, x, par(1),par(2),par(3),par(4),par(5),par(6),par(7),par(8),par(9),par(10),par(11),par(12),par(13),par(14));
    [xs,~,flag] = fsolve(F, seed, opts_fs);
    if flag <= 0, out = NaN; return; end
    J = jacobian(0, xs, par(1),par(2),par(3),par(4),par(5),par(6),par(7),par(8),par(9),par(10),par(11),par(12),par(13),par(14));
    out = sign(trace(J));
end

%%% --------------------------------- Jacobians and hessians ---------------------------------
function jac = jacobian(~,x,w_EE,w_EI,w_IE,w_II,P_E,P_I,k_E,x0_E,k_I,x0_I,rE,rI,tau_E,tau_I)
    E=x(1); I=x(2);
    xE = w_EE*E - w_EI*I + P_E;
    xI = w_IE*E - w_II*I + P_I;
    SE = 1/(1+exp(-k_E*(xE-x0_E)));
    SI = 1/(1+exp(-k_I*(xI-x0_I)));
    cE = 1/(1+exp(k_E*x0_E));
    cI = 1/(1+exp(k_I*x0_I));
    dSE = k_E*SE*(1-SE);
    dSI = k_I*SI*(1-SI);
    dfE_dE = (-1 - rE*(SE-cE) + (1-rE*E)*dSE*w_EE)/tau_E;
    dfE_dI = (-(1-rE*E)*dSE*w_EI)/tau_E;
    dfI_dE = ((1-rI*I)*dSI*w_IE)/tau_I;
    dfI_dI = (-1 - rI*(SI-cI) - (1-rI*I)*dSI*w_II)/tau_I;
    jac = [dfE_dE dfE_dI; dfI_dE dfI_dI];
end

function jacp = jacobianp(~,x,w_EE,w_EI,w_IE,w_II,P_E,P_I,k_E,x0_E,k_I,x0_I,rE,rI,tau_E,tau_I)
    E=x(1); I=x(2);
    xE = w_EE*E - w_EI*I + P_E;
    xI = w_IE*E - w_II*I + P_I;
    SE = 1/(1+exp(-k_E*(xE-x0_E)));
    SI = 1/(1+exp(-k_I*(xI-x0_I)));
    cE = 1/(1+exp(k_E*x0_E));
    cI = 1/(1+exp(k_I*x0_I));
    fE = (-E + (1-rE*E)*(SE-cE))/tau_E;
    fI = (-I + (1-rI*I)*(SI-cI))/tau_I;

    dfE_dwEE = (1-rE*E)*SE*(1-SE)*k_E*E/tau_E;
    dfE_dwEI = -k_E*I*(1-rE*E)*SE*(1-SE)/tau_E;
    dfE_dPE  = (1-rE*E)*SE*(1-SE)*k_E/tau_E;
    dfE_dkE  = (1-rE*E)*( SE*(1-SE)*(xE-x0_E) + cE*(1-cE)*x0_E )/tau_E;
    dfE_dx0E = (1-rE*E)*k_E*( cE*(1-cE) - SE*(1-SE) )/tau_E;
    dfE_drE  = (-E*(SE-cE))/tau_E;
    dfE_dtauE= -fE/tau_E;

    dfI_dwIE = (1-rI*I)*SI*(1-SI)*k_I*E/tau_I;
    dfI_dwII = -k_I*I*(1-rI*I)*SI*(1-SI)/tau_I;
    dfI_dPI  = (1-rI*I)*SI*(1-SI)*k_I/tau_I;
    dfI_dkI  = (1-rI*I)*( SI*(1-SI)*(xI-x0_I) + cI*(1-cI)*x0_I )/tau_I;
    dfI_dx0I = (1-rI*I)*k_I*( cI*(1-cI) - SI*(1-SI) )/tau_I;
    dfI_drI  = (-I*(SI-cI))/tau_I;
    dfI_dtauI= -fI/tau_I;

    jacp = [ dfE_dwEE dfE_dwEI 0 0 dfE_dPE 0 dfE_dkE dfE_dx0E 0 0 dfE_drE 0 dfE_dtauE 0 ;
             0 0 dfI_dwIE dfI_dwII 0 dfI_dPI 0 0 dfI_dkI dfI_dx0I 0 dfI_drI 0 dfI_dtauI ];
end

function hess = hessians(~,x,w_EE,w_EI,w_IE,w_II,P_E,P_I,k_E,x0_E,k_I,x0_I,rE,rI,tau_E,tau_I)
    E=x(1); I=x(2);
    xE = w_EE*E - w_EI*I + P_E;
    xI = w_IE*E - w_II*I + P_I;
    SE = 1/(1+exp(-k_E*(xE-x0_E)));
    SI = 1/(1+exp(-k_I*(xI-x0_I)));
    dSE = k_E*SE*(1-SE);
    dSI = k_I*SI*(1-SI);
    d2SE = k_E^2*SE*(1-SE)*(1-2*SE);
    d2SI = k_I^2*SI*(1-SI)*(1-2*SI);

    d2fE_dEdE = ( -2*rE*w_EE*dSE + (1-rE*E)*w_EE^2*d2SE )/tau_E;
    d2fE_dEdI = (  rE*w_EI*dSE  - (1-rE*E)*w_EE*w_EI*d2SE )/tau_E;
    d2fE_dIdI = ( (1-rE*E)*w_EI^2*d2SE )/tau_E;

    d2fI_dEdE = ( (1-rI*I)*w_IE^2*d2SI )/tau_I;
    d2fI_dEdI = ( -rI*w_IE*dSI - (1-rI*I)*w_IE*w_II*d2SI )/tau_I;
    d2fI_dIdI = (  2*rI*w_II*dSI + (1-rI*I)*w_II^2*d2SI )/tau_I;

    hess(:,:,1) = [d2fE_dEdE d2fE_dEdI; d2fE_dEdI d2fE_dIdI];
    hess(:,:,2) = [d2fI_dEdE d2fI_dEdI; d2fI_dEdI d2fI_dIdI];
end
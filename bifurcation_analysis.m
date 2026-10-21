%% Bifurcation analysis of wilson Cowan model
%
% System:
%   dE/dt = (-E + (1 - p.rE*E) / (1 + exp(-p.k_E * (xE - p.x0_E))) - 1/(1 + exp(p.k_E * p.x0_E)))/p.tau_E
%   dI/dt = (-I + (1 - p.rI*I) / (1 + exp(-p.k_I * (xI - p.x0_I))) - 1/(1 + exp(p.k_I * p.x0_I)))/p.tau_I

%% Step 1: Define parameters and initial equilibrium
clear;
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
cpl(x,v,s,[3,1]);

x1 = x(1:2, s(2).index);
par(ap1) = x(end, s(2).index);
[x0,v0] = init_H_H(@continuation_equations, x1, par, [2,4]);

opt = contset(opt, 'Backward', 0);
[x2f,v2f,s2f,h2f,f2f] = cont(@hopf, x0, v0, opt);

opt = contset(opt, 'Backward', 1);
[x2b,v2b,s2b,h2b,f2b] = cont(@hopf, x0, v0, opt);

figure; hold on;
plot(x2f(3,:), x2f(4,:), 'b-', 'LineWidth', 1.8);
plot(x2b(3,:), x2b(4,:), 'b-', 'LineWidth', 1.8);
plot(x0(3), x0(4), 'ko', 'MarkerFaceColor','r', 'MarkerSize',7);  % starting H point
xlabel('w_{EI}', 'FontSize', 12);
ylabel('w_{IE}', 'FontSize', 12);
title('Hopf curve in (w_{EI}, w_{IE})', 'FontSize', 13);

k = 2;                               % the index where label is 'H'
x1 = x(1:2, s(k).index);
par(ap1) = x(end, s(k).index);
[x0,v0] = init_H_H(@continuation_equations, x1, par, [2,4]);   % w_EI, w_II

opt = contset(opt, 'Backward', 0);
[x2f,v2f,s2f] = cont(@hopf, x0, v0, opt);
opt = contset(opt, 'Backward', 1);
[x2b,v2b,s2b] = cont(@hopf, x0, v0, opt);

figure; hold on;
plot(x2f(3,:), x2f(4,:), 'b-', 'LineWidth', 1.8);
plot(x2b(3,:), x2b(4,:), 'b-', 'LineWidth', 1.8);
plot(x0(3), x0(4), 'ko', 'MarkerFaceColor','r', 'MarkerSize',7);
xlabel('w_{EI}');
ylabel('w_{II}');
title('Hopf curve in (w_{EI}, w_{II}) — decreasing branch');
grid on; box on;

% combine both branches of the curve into one ordered list
wEI_c = [fliplr(x2b(3,:)) x2f(3,:)];
wII_c = [fliplr(x2b(4,:)) x2f(4,:)];
Eq_c  = [fliplr(x2b(1:2,:)') ; x2f(1:2,:)'];

delta = 0.6;  % small perpendicular offset in w_II
side_in  = nan(1,length(wEI_c));
side_out = nan(1,length(wEI_c));
for i = 1:length(wEI_c)
    seed = Eq_c(i,:)';
    side_in(i)  = side_class(wEI_c(i), wII_c(i)-delta, seed, par, opts_fs);
    side_out(i) = side_class(wEI_c(i), wII_c(i)+delta, seed, par, opts_fs);
end

figure; hold on;
plot(wEI_c(side_in>0),  wII_c(side_in>0)-delta,  '.', 'Color',[1 0.8 0.75], 'MarkerSize',10);
plot(wEI_c(side_in>0),  wII_c(side_in>0)-delta,  '.', 'Color',[1 0.8 0.75], 'MarkerSize',10);
plot(wEI_c, wII_c, 'k-', 'LineWidth', 2);
plot(x0(3), x0(4), 'ko', 'MarkerFaceColor','r', 'MarkerSize',7);
xlabel('w_{EI}'); ylabel('w_{II}');
title('Hopf curve in (w_{EI}, w_{II})');
grid on; box on;


wEI_c = [fliplr(x2b(3,:)) x2f(3,:)];
wII_c = [fliplr(x2b(4,:)) x2f(4,:)];
Eq_c  = [fliplr(x2b(1:2,:)') ; x2f(1:2,:)'];

delta = 0.6;
below_pts = []; above_pts = [];
for i = 1:length(wEI_c)
    seed = Eq_c(i,:)';
    s_below = side_class(wEI_c(i), wII_c(i)-delta, seed, par, opts_fs);
    s_above = side_class(wEI_c(i), wII_c(i)+delta, seed, par, opts_fs);
    if ~isnan(s_below), below_pts = [below_pts; wEI_c(i), wII_c(i)-delta, s_below]; end
    if ~isnan(s_above), above_pts = [above_pts; wEI_c(i), wII_c(i)+delta, s_above]; end
end

fprintf('below curve: trace sign = %+.0f (consistent in %d%% of samples)\n', mode(below_pts(:,3)), round(100*mean(below_pts(:,3)==mode(below_pts(:,3)))));
fprintf('above curve: trace sign = %+.0f (consistent in %d%% of samples)\n', mode(above_pts(:,3)), round(100*mean(above_pts(:,3)==mode(above_pts(:,3)))));
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
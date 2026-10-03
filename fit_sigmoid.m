function [k, x0, r_ref] = fit_sigmoid(I_list, nu, I_max,type)
    % FIT_SIGMOID  Original Wilson-Cowan sigmoid  F(x) = 1 / (1 + exp(-k (x - x0)))
    %
    %   [k, x0, r_ref] = fit_sigmoid(I_list, nu, I_max,type)
    %
    % Inputs:
    %   I_list : injected currents (uA/cm^2)
    %   nu     : firing rates (Hz), same length as I_list
    %   I_max  : upper end of fit range (uA/cm^2)   [optional, default: all]
    %
    % Outputs:
    %   k      : slope
    %   x0     : threshold (uA/cm^2)
    %   r_ref  : fixed rate scale (Hz) — rates are fitted as nu / r_ref
    
    I_list = I_list(:);  nu = nu(:);
    if nargin < 3, I_max = max(I_list); end
    
    % ---------- fit range ----------
    keep = I_list >= 0 & I_list <= I_max;
    x = I_list(keep);
    
    % ---------- normalize rates to 0..1 ----------
    r_ref = max(nu(keep));
    y     = nu(keep) / r_ref;
    
    % ---------- fit k and x0 ----------
    sig = @(q, x) 1 ./ (1 + exp(-abs(q(1)) * (x - q(2))));      % q = [k, x0]
    q0  = [4/(max(x) - min(x)), x(find(y >= 0.5, 1))];
    q   = fminsearch(@(q) sum((sig(q, x) - y).^2), q0, optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5, 'TolFun', 1e-12, 'TolX', 1e-12));
    
    k  = abs(q(1));
    x0 = q(2);
    
    % ---------- plot ----------
    figure;
    plot(x, y, 'o', x, sig([k x0], x), 'k-', 'LineWidth', 1.5);
    xlabel('I [\muA/cm^2]'); ylabel('rate / r_{ref}');
    legend('neuron', 'sigmoid fit', 'Location', 'northwest');
    title(sprintf('k = %.3f,  x_0 = %.3f,  r_{ref} = %.1f Hz, type = %s', k, x0, r_ref,type));
    fprintf('For: %s, k = %.3f,  x_0 = %.3f, r_ref = %.3f \n', type, k, x0,r_ref)
end
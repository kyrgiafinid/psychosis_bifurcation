function dy = wilson_cowan(t, y, p)
    % WILSON_COWAN  Right-hand side of the Wilson-Cowan (1972) model
    %
    %   tau_E dE/dt = -E + (1 - rE*E) * S_E( w_EE*E - w_EI*I + P_E )
    %   tau_I dI/dt = -I + (1 - rI*I) * S_I( w_IE*E - w_II*I + P_I )
    %
    %   S(x) = 1 / (1 + exp(-k (x - x0)))
    %
    % y = [E; I]   activity of the PY (E) and IN (I) populations
    % p = struct with fields
    %     w_EE, w_EI, w_IE, w_II   weights
    %     P_E, P_I                 external input
    %     k_E, x0_E, k_I, x0_I     sigmoid slope and threshold
    %     rE, rI                   refractory terms
    %     tau_E, tau_I             time constants (ms)
    
    E = y(1);
    I = y(2);
    
    % --- total input to each population ---
    xE = p.w_EE*E - p.w_EI*I + p.P_E;
    xI = p.w_IE*E - p.w_II*I + p.P_I;
    
    % --- sigmoid response ---
    SE = 1 / (1 + exp(-p.k_E * (xE - p.x0_E)));
    SI = 1 / (1 + exp(-p.k_I * (xI - p.x0_I)));
    
    % --- fraction of cells not refractory ---
    avail_E = max(0, 1 - p.rE*E);
    avail_I = max(0, 1 - p.rI*I);
    
    % --- equations ---
    dy = [(-E + avail_E * SE) / p.tau_E;
          (-I + avail_I * SI) / p.tau_I];

end
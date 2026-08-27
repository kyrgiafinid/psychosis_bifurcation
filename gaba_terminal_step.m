function [X,Y,c,b] = gaba_terminal_step(X,Y,c,b,spikedIN,dt,U,Util,tau_D,tau_R,lambda_t,Ka,alpha_ca,KP,beta_ca,kplus,kminus,bt,IP)
    Z = 1 - X - Y;
    X = X + dt.*(Z./tau_R);
    Y = Y + dt.*(-Y./tau_D);

    lam = lambda_t .* (c.^4) ./ (c.^4 + Ka.^4);
    asyncEvent = rand(size(c)) < (1 - exp(-lam.*dt));

    pump = alpha_ca .* c.^2 ./ (c.^2 + KP.^2);
    buffer = kplus.*(bt-b).*c - kminus.*b;
    c = c + dt.*(-pump - buffer + IP);
    b = b + dt.*buffer;

    if any(spikedIN)
        used = U.*X(spikedIN);
        X(spikedIN) = X(spikedIN) - used;
        Y(spikedIN) = Y(spikedIN) + used;
        c(spikedIN) = c(spikedIN) + beta_ca;
    end
    if any(asyncEvent)
        used = Util.*X(asyncEvent);
        X(asyncEvent) = X(asyncEvent) - used;
        Y(asyncEvent) = Y(asyncEvent) + used;
    end
    X = min(max(X,0),1); Y = min(max(Y,0),1);
    c = max(c,1e-6); b = min(max(b,0),bt);
end
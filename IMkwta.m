function [B, committed] = IMkwta(B, duration)
% implements k-WTA through global inhibition with the Runge-Kutta 4 method

global P
global C

GI = @(x) P.selfActB.*(P.asyB-x) - P.inhibB*sum(x(:)); % non-shunting version
nSteps = round(duration./C.tstep);

for t = 1:nSteps
    % RK4
    k1 = GI(B);
    k2 = GI(max(0, min(P.asyB, B + 0.5*C.tstep*k1)));
    k3 = GI(max(0, min(P.asyB, B + 0.5*C.tstep*k2)));
    k4 = GI(max(0, min(P.asyB, B + C.tstep*k3)));
    B = max(0, min(P.asyB, B + (C.tstep/6) * (k1 + 2*k2 + 2*k3 + k4)));
end
committed = B > 0;


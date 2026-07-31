function Committed = IMkwtaPre(maxDuration)
% implements k-WTA through global inhibition with the Runge-Kutta 4 method
% used to pre-calculate all combinations of initial recruitment and release
% duration

global P
global C

GI = @(x) P.selfActB.*(P.asyB-x) - P.inhibB*sum(x(:)); % non-shunting version
Timeline = C.tstep:C.tstep:maxDuration;
nRuns = 30;
Committed = zeros(P.nb, length(Timeline));
tic
for run = 1:nRuns
for initRecruit = 1:P.nb
    B = exp(randn(1, initRecruit))./sqrt(200); % Initialize B for each initial recruitment
    for duration = Timeline
        nSteps = round(duration./C.tstep);
        for t = 1:nSteps
            % RK4
            k1 = GI(B);
            k2 = GI(max(0, min(P.asyB, B + 0.5*C.tstep*k1)));
            k3 = GI(max(0, min(P.asyB, B + 0.5*C.tstep*k2)));
            k4 = GI(max(0, min(P.asyB, B + C.tstep*k3)));
            B = max(0, min(P.asyB, B + (C.tstep/6) * (k1 + 2*k2 + 2*k3 + k4)));
        end
        Committed(initRecruit, nSteps) = Committed(initRecruit, nSteps) + sum(B > 0);
    end
end
end
Committed = Committed./nRuns;
toc

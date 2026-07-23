function [Map, W] = IMdecayFX(Map, W, delay, tstep, input)
% decay of the feature map and the weight matrix
global C
global P

if nargin < 4, tstep = C.tstep; end
if nargin < 5, input = zeros(size(Map(1).FX)); end

Decay = @(x) P.selfactFX.*x - P.inhibFX*sum(x(:)) + P.stimDrive*(P.asyFX-max(x(:)))*input;

timesteps = round(max(0,delay)/tstep);  % delay is in s
for ff = 1:C.nfeatures
    for t = 1:timesteps  

        %summedAct = sum(Map(ff).FX(:));
        %Map(ff).FX = max(0, min(P.asyFX, P.selfactFX*Map(ff).FX) - P.inhibFX*summedAct);

        FX = Map(ff).FX;
        k1 = Decay(FX);
        k2 = Decay(max(0, min(P.asyFX, FX + 0.5*tstep*k1)));
        k3 = Decay(max(0, min(P.asyFX, FX + 0.5*tstep*k2)));
        k4 = Decay(max(0, min(P.asyFX, FX + tstep*k3)));
        Map(ff).FX = max(0, min(P.asyFX, FX + (tstep/6) * (k1 + 2*k2 + 2*k3 + k4)));

        % Map(ff).FX = DecayFX_mex(Map(ff).FX, delay, P.selfactFX, P.inhibFX, P.asyFX, tstep);
        % The mex is actually slower!

    end
end

W = W + randn(size(W)) * sqrt(max(0, delay)*P.wnoise.^2); % variance is added at a constant rate over time --> sqrt takes the SD
output = 1; 



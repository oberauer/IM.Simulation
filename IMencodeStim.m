function [W, GateClosed, GateWeight, committedNew, bStrength] = IMencodeStim(W, context, content, GateClosed, GateWeight, consRate, consTime, releaseTime)

global P
global C

contentTransposed = content';
inputVec = [context, contentTransposed(:)'];   % the vectorization of the now vertical contents concatenates all features of multi-feature items
freeBP = rand(1, P.nb) < (1-P.delta);
GateClosed(freeBP) = 0;
W(:, freeBP) = 0;  % remove weights to the now free binding units

strength = 1 - exp(-consRate*consTime); 
%strength = 1; 

initRecruited = rand(1, P.nb) < P.pMax;
Binding = randn(1, P.nb)./P.nbNorm;   % normalize by the average norm of nb across individuals, based on mean P.nb across individuals (and conditions, when P.nb is experimentally varied)
Binding = Binding .* (1-GateClosed) .* initRecruited; 
[Bsorted, sortIdx] = sort(abs(Binding), 'descend'); 
ncommit = round( sum((1-GateClosed) .* initRecruited) * (P.pBase + (1-P.pBase) * exp(-P.rRate*releaseTime) ) );  % number of binding units that will remain committed -> gate will be closed
committedNew = sortIdx(1:ncommit);
BindingsNow = zeros(1, P.nb);
BindingsNow(committedNew) = Binding(committedNew); 
GateClosed(committedNew) = 1;  % take the nbind free units with the highest absolute gating values and close the gate for them (= commit them)
GateWeight(committedNew) = GateWeight(committedNew) + abs(Binding(committedNew));
bStrength = abs(BindingsNow); 
W = W + strength * (inputVec' * BindingsNow); 
output = W; 


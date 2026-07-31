function [W, GateClosed, GateWeight, committedNew, bStrength] = IMencodeStim(W, context, content, GateClosed, GateWeight, consRate, consTime, releaseTime)

global P

contentTransposed = content';
inputVec = [context, contentTransposed(:)'];   % the vectorization of the now vertical contents concatenates all features of multi-feature items
freeBP = rand(1, P.nb) < (1-P.delta);
GateClosed(freeBP) = 0;
W(:, freeBP) = 0;  % remove weights to the now free binding units
strength = 1 - exp(-consRate*consTime); 
Binding = exp(randn(1, P.nb))./P.nbNorm;
Binding = Binding .* (1-GateClosed);   % set the committed units to zero 

nRecruited = sum(GateClosed);
ncommit = C.KWTA(nRecruited, round(min(releaseTime,1)/C.tstep));
[Bsorted, sortIdx] = sort(abs(Binding), 'descend'); 
committedNew = sortIdx(1:ncommit);
BindingsNow = zeros(1, P.nb);
BindingsNow(committedNew) = Binding(committedNew); 
GateClosed(committedNew) = 1;  % take the newly committed binding units and close the gate for them (= commit them)
GateWeight(committedNew) = GateWeight(committedNew) + abs(Binding(committedNew));
bStrength = abs(BindingsNow); 
W = W + strength * (inputVec' * BindingsNow); 
output = W; 

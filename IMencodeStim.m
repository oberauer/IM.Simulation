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
[Binding, committedNew] = IMkwta(Binding, releaseTime); 
GateClosed(committedNew) = 1;  
GateWeight(committedNew) = GateWeight(committedNew) + abs(Binding(committedNew));
bStrength = Binding; 
W = W + strength * (inputVec' * Binding); 
output = W; 


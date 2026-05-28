function W = CreateConnections(nfeat)

global C
global P
global E

if nargin==0, nfeat=1; end

% mapping from categories to binding pool

if E.context == 1, nCtxCat = C.nLocCat; end
if E.context == 2, nCtxCat = C.nCat; end

% connection weight matrix for all inputs (context and all the features, concatenated) to the binding layer

if E.context == 1, W = zeros(nCtxCat + nfeat*C.nCat, P.nb); end
if E.context == 2, W = zeros(nfeat*C.nCat, P.nb); end   % only the cue feature and the target feature are bound together, the location is ignored. 

wait = 1; 



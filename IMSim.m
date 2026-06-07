function output = IMSim(P, setsize, cueing, refreshings, MemoryState)
% Function for simulating one trial of the IM model of visual WM.

global E
global C

%%%%%%% Preparatory steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    if cueing == 4, refreshings = 1; else, refreshings = 0; end % by default, assume no refreshing during retention interval 
end
if nargin < 5
    map = struct('FX', zeros(C.nc));   % feature map
    Map = repmat(map, C.nfeatures, 1);
    W = CreateConnections(C.nfeatures);
    G = zeros(1, P.nb);  % gate-closing units
    GW = zeros(1, P.nb); % gate-closing weights
else
    Map = MemoryState{1};
    W = MemoryState{2};
    G = MemoryState{3}; 
    GW = zeros(1, P.nb); % gate-closing weights
end

% Generate memory set for this trial
if E.layout == 1, L = randperm(C.nloc); end      %shuffle locations of array objects
if E.layout == 2, L = ones(1, C.nloc); end       % all in location 1
if E.layout == 3, L = randperm(floor(C.nloc/2)); end  % all on the first semi-circle
if E.layout == 4, L = randperm(floor(C.nloc/2)) + ceil(C.nloc/2); end  % all on the second semi-circle
F = zeros(C.nfeatures, C.nstim);
for ff = 1:C.nfeatures
    F(ff,:) = randperm(C.nstim);      %shuffle object features
end

if ismember(E.test, 2:4)  % if this is a change-detection (=recognition) or n-AFC or change localization test
    [probestim, probeIdx] = IMprepareProbe(F, setsize);
else
    probestim = []; probeIdx = []; 
end

%%%%%%%%%%% Encoding and Consolidation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Map, W, G, GW, Focus, Afocus, content, context, Inpos, Strength, Bstrength, CTime, SpatAttn] = IMencoding(Map, W, G, GW, L, F, setsize, cueing);
usedTime = setsize*CTime + P.maskWindow; %CTime is the mean consolidation time taken
overTime = max(0, usedTime - (E.RI + E.prestime));

g1 = GW;  % keep record of "gate closed" BP units (GW: with their G weight)
fx1 = Map(1).FX;

% if E.presentation == 1, delay = E.prestime - sum(cTime) + E.RI; end
% if E.presentation == 2, delay = E.prestime - cTime(setsize) + E.RI; end
%Map = IMdecayFX(Map, delay);  % decay of FX
%w = w - (1-P.delta)*delay*w;  % decay of w

%%%%%%%%%%%% Now do all the retro-cueing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if cueing > 1 && E.PreRetro == 2
    [Map, W, G, GW, Focus, Afocus, lastrefreshed, CuedIdx, Strength] = IMcueing(Map, W, G, GW, Strength, Focus, Afocus, L, F, setsize, cueing, refreshings);
    %disp(max(abs(W(:))));
else
    lastrefreshed = 0;
end

str = Strength(1:E.outsize); % content strength of the target
g = G;    % record the state of the GateClosed vector for carry-over to a second to-be-encoded array
gw = GW;  % record the state of GateWeights before test (during which it will be changed by output interference)
w = W;  % record the state of W before test (during which it will be changed by output interference and/or decay)
fx = Map(1).FX;

%%%%%%%%%%%% Test: recall or recognition ##########################

if E.test == 1 || E.test == 4, response = zeros(1, E.outsize); end
if E.test == 2 || E.test == 3, response = zeros(2, E.outsize); end
rt = zeros(1,E.outsize);
for probed = 1:E.outsize
    [response(:,probed), rt(probed), Map, W, G, Focus, CWcolor] = IMretrieve(Map, W, G, Focus, Afocus, probed, cueing, L, F, probestim, probeIdx, overTime);
end  % outpos

%%%%%%%%%%%%%%% Wrap up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% put together everything we want the function to return in a structure
output.response = response;
output.rt = rt;
output.F = F;
output.L = L;
output.content = content;
output.context = context; 
output.lastrefreshed = lastrefreshed;
output.wx = w; 
output.g = g;
output.fx = fx;
output.map = Map;
output.Inpos = Inpos;
output.Strength = str;
output.Bstrength = Bstrength; 
output.CTime = CTime;
output.SpatAttn = SpatAttn;
output.CWcolor = CWcolor;

if E.test == 2   % additional outputs for recognition
    output.SameRange = E.SameRange;
    output.probeIdx = probeIdx;   %location of origin of intrusion probes
    if cueing == 6, output.cueIdx = CuedIdx; end
end


end
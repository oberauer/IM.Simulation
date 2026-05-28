function output = IMSimAlpha(P, setsize, plotTrajectory)
% Function for simulating one trial of the IM model of visual WM, recording
% the time course of spatial attention as modulated by alpha wave

global E
global C

map = struct('FX', zeros(C.nc));   % feature map
Map = repmat(map, C.nfeatures, 1);
B = zeros(1, P.nb);  % binding pool
G = zeros(1, P.nb);  % gate-closing units
GW = zeros(1, P.nb); % gate-closing weights

% Generate memory set for this trial
if E.layout == 1, L = randperm(C.nloc, E.maxsetsize); end %shuffle locations of array objects
if E.layout == 2, L = ones(1, E.maxsetsize); end  % all in location 1
F = zeros(C.nfeatures, C.nstim);
for ff = 1:C.nfeatures
    F(ff,:) = randperm(C.nstim);      %shuffle object features
end

if E.test == 2  % if this is a change-detection (=recognition) test
    [probestim, probeIdx, SameRange, ChangeRange] = IMprepareRecog(F, setsize);
else
    probestim = []; SameRange = []; ChangeRange = []; probeIdx = [];
end


%%%%%%%%%%% Encoding and Consolidation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Map, B, G, GW, Focus, Afocus, content, context, Inpos, Strength, Bstrength, SpatAttn] = IMencodingAlpha(Map, B, G, GW, L, F, setsize, plotTrajectory);

b = B;
g = GW;  % keep record of "gate closed" BP units (GW: with their G weight)
fx = Map(1).FX;

%%%%%%%%%%%% Test: recall or recognition ##########################

[Testdisplay, CWcolor] = IMdisplay(L, probestim);

response = zeros(E.test,E.outsize);
rt = zeros(1,E.outsize);
for probed = 1:E.outsize
    [response(:,probed), rt(probed), Map, B, G, Focus] = IMretrieve(Map, B, G, Focus, Afocus, probed, 1, L, F, Testdisplay, SameRange, ChangeRange, probeIdx);
end  % outpos

%%%%%%%%%%%%%%% Wrap up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% put together everything we want the function to return in a structure
output.response = response;
output.rt = rt;
output.F = F;
output.L = L;
output.content = content;
output.context = context; 
output.lastrefreshed = 0;
output.b = b; 
output.g = g;
output.fx = fx;
output.map = Map;
output.Inpos = Inpos;
output.Strength = Strength;
output.Bstrength = Bstrength;
output.SpatAttn = SpatAttn;
output.CWcolor = CWcolor;

if E.test == 2   % additional outputs for recognition
    output.SameRange = SameRange;
    output.probeIdx = probeIdx;   %location of origin of intrusion probes
    if cueing == 6, output.cueIdx = CuedIdx; end
end


end
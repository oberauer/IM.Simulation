function [Map, W, CDAg, CDAw, Alpha, sumCtime, Pangle, EEG_W, EEG_FX] = IMtrackSignals(setsize, eW, eW2, eNoise)
% encodes a memory set simultaneously, simulates the proces time-step wise
% to track CDA and alpha over time

global P
global E
global C

map = struct('FX', zeros(C.nc));   % feature map
Map = repmat(map, C.nfeatures, 1);
W = CreateConnections(C.nfeatures);
G = zeros(1, P.nb);
GW = zeros(1, P.nb);

% Generate memory set for this trial
if E.layout == 1, L = randperm(C.nloc); end      %shuffle locations of array objects
if E.layout == 2, L = ones(1, C.nloc); end       % all in location 1
F = zeros(C.nfeatures, C.nstim);
for ff = 1:C.nfeatures
    F(ff,:) = randperm(C.nstim);      %shuffle object features
end

if E.forwardrecall, encorder = 1:setsize; else, encorder = randperm(setsize); end % for consolidation: array items are encoded sequentially (in a random order chosen by the subject). This variable codes the planned output position of the encoded items

% determine whether a stimulus is integrated with a mask
maskWindow = ones(1, setsize) * P.maskWindow;
sdmaskWindow = P.maskWindowSD * P.maskWindow;
attentionWindows = gamrnd(maskWindow.^2./sdmaskWindow.^2, sdmaskWindow.^2./maskWindow, 1, setsize);
masking = attentionWindows > abs(E.MaskSOA);

cRate = gamrnd(P.cRate^2/P.cRateSD^2, P.cRateSD^2/P.cRate, 1, setsize+1);  % generate the vector of consolidation rates
cTime = -(log(1-P.cStrength)./cRate); % vector of consolidation times needed to reach strength P.cTau
cTime = diff([0, min(E.RI, attentionWindows(encorder(1)) + cumsum(cTime))]);

Timepoints = E.RI/C.tstep + 1;
CDAg = zeros(1,Timepoints);
CDAw = zeros(1,Timepoints);
Alpha = zeros(1,Timepoints);
nElectrodes = size(eW, 2);
EEG_W = zeros(Timepoints, nElectrodes);
EEG_FX = zeros(Timepoints, nElectrodes);

%sequential encoding - binding (read out from FX in sequential order)
Inpos = zeros(1,setsize);  % initialize vector coding the input positions of the items (in the encoding order).
% In this vector, and the result vectors, the items are ordered by output position!

%Focus = randperm(setsize, 1);  % start with random location of the focus
% AfocusLoc = C.location(L(Focus),:);     % location attended to
% Afocus = zeros(1,C.nc);                 % just to initialize it
SpatAttn = zeros(C.nc, 1);

% simultaneous presentation: parallel encoding into spatially organized feature maps;
% replacment by mask if the mask falls within the replacement window
for ff = 1:C.nfeatures
    for inpos = 1:setsize
        Map(ff).FX = Map(ff).FX + masking(inpos) .* C.maskStim + (1-masking(inpos)) .* C.location(L(inpos),:)' * (C.stim(F(ff,inpos),:));
    end
end

cumCtime = [0, cumsum(cTime(1:setsize)), E.RI];
sumCtime = sum(cTime);
inpos = 0;
t = 0;
tcount = 1;
consolStarted = zeros(1, setsize);

while t < E.RI  % continue until end of RI
    % presentation of next stimulus: encode into feature map
    if t >= cumCtime(inpos+1)  % if the cTime for item inpos+1 starts, ...
        inpos = inpos + 1;  % ... move to inpos+1
    end

    if inpos <= setsize
        if consolStarted(inpos) == 0
            % initial consolidation (for 1 time step)
            Focus = encorder(inpos);
            AfocusLoc = C.location(L(Focus),:);     % update location attended to in the feature maps
            Afocus = AfocusLoc./sum(AfocusLoc) * Map(1).FX; % use location as (spatial) attentional filter to pull out the target feature from its feature map
            content = C.stimnoise + Afocus * C.Mapping;
            if E.context == 1, context = C.locationnoise + AfocusLoc * C.MappingC; end
            if E.context == 2, context = C.stimnoise + (AfocusLoc./sum(AfocusLoc) * Map(2).FX) * C.Mapping; end
            [W, G, GW, committedNew, ~] = IMencodeStim(W, context, content, G, GW, cRate(inpos), C.tstep, C.tstep);
            committedNewNotBase = randsample(committedNew, round(length(committedNew)*(1-P.pBase))); % sample a subset of newly committed binding units as the to-be-released ones
            consolStarted(inpos) = 1;
        else
            % continuing consolidation
            pLoss = 1 - 1/exp(cRate(inpos)./(cTime(inpos).*C.tstep)); % see StepByStepEncoding.m
            nLoss = binornd(length(committedNewNotBase), pLoss);
            decommitted = randsample(committedNewNotBase, nLoss);
            committedNewNotBase = setdiff(committedNewNotBase, decommitted);
            G(decommitted) = 0;
            W(:, decommitted) = 0;  % remove weights to the now free binding units
        end
    end

    [Map, W] = IMdecayFX(Map, W, C.tstep);   % decay of FX through one time step
    SpatAttn = max( 0, SpatAttn + C.tstep * (mean(Map(1).FX,2) + P.TopDownSpatAttn.*AfocusLoc' - P.spatinhib*SpatAttn*sum(SpatAttn)) );  % attraction of spatial attention to locations in feature maps, top-down guidance by FoA, and global inhibition on spatial attention
    %Map(1).FX = Map(1).FX + C.tstep * (-Map(1).FX + Map(1).FX .* repmat(SpatAttn, 1, C.nc)); % spatial attention modulates feature maps

    EEG_W(tcount,:) = (((context * W(1:C.nLocCat, :)) * W((C.nLocCat+1):end, :)') * C.Mapping') * eW + randn(1,nElectrodes)*eNoise;  % feed last-used context into weight matrix -> reactivate content -> project onto electrodes
    EEG_FX(tcount,:) = SpatAttn' * eW + randn(1,nElectrodes)*eNoise;
    CDAg(tcount) = sum(G);
    CDAw(tcount) = sum(abs(W(:)));
    Alpha(tcount) = sum(abs(EEG_FX(tcount,:)));
    tcount = tcount + 1;
    t = t + C.tstep;
end

for inpos = 1:setsize
    Inpos(encorder(inpos)) = inpos;
end
Pangle = L(1:setsize);

end


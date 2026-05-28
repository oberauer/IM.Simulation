function [Map, B, GateClosed, gWeight, Focus, Afocus, content, context, Inpos, Strength, Bstrength, SpatAttnRecord] = IMencodingAlpha(Map, B, GateClosed, gWeight, L, F, setsize, plotTrajectory)
% encodes a simultaneously presented memory array by circulating the FoA
% among items in the FX map in an alpha rhythm

global P
global E
global C

if nargin < 8, plotTrajectory = 0; end
P.IOR = 0.02; % strength of inhibition of return

P.spatinhib = 4;
minInhib = P.spatinhib*0.05; 
tconst = 0.05; % time constant of global-inhibition unit in s
tlag = round(tconst/C.tstep); % time constant in time steps
attnBaseInput = 0.2; 
attnThreshold = attnBaseInput + 0.15;
maxAttn = 0.02;

cRate = gamrnd(P.cRate^2/P.cRateSD^2, P.cRateSD^2/P.cRate, 1, 100);  % generate the vector of consolidation rates
cTime = -(log(P.pBase)./cRate); % vector of consolidation times needed to reach binding-pool recruitment pBase
        
Inpos = zeros(1,setsize);  % initialize vector coding the input positions of the items (in the encoding order). In this vector, and the result vectors, the items are ordered by output position!
Focus = randperm(setsize, 1);  % start with random location of the focus
AfocusLoc = zeros(1,C.nc);
Afocus = zeros(1,C.nc);                 % just to initialize it
committed = zeros(setsize, P.nb);   % keeps a record of which binding units are committed (i.e., gate closed) for each item
bstrength = zeros(setsize, P.nb);

% simultaneous presentation: parallel encoding into spatially organized feature maps
for ff = 1:C.nfeatures
    for inpos = 1:setsize
        Map(ff).FX = Map(ff).FX + C.location(L(inpos),:)' * C.stim(F(ff,inpos),:);
    end
end
%SpatAttn = (mean(Map(1).FX,2)) + randn(C.nc, 1) * 0.0001;        % start with some random noise
SpatAttn = zeros(C.nc, 1);
SpatAttnRecord = zeros(ceil(E.RI/C.tstep), C.nc);
SpatInhib = 0;
SpatInhibRecord = zeros(1, E.RI/C.tstep); 

cClock = 0; % clock for current consolidation process
t = C.tstep;
tIdx = 1;
inpos = 1;
cNum = 1;
currentCons = 0; % which item is currently being consolidated
consSequence = zeros(1,20);
startNewC = 1;
Strength = zeros(1,setsize);
FStrength = zeros(E.RI/C.tstep, setsize); % for tracking strength of features in FX
SStrength = zeros(E.RI/C.tstep, setsize); % for tracking strength of spatial attention foci in FX
tx = 1;
inhibVec = P.spatinhib*ones(C.nc, 1); 

while t < E.RI  % continue until end of presentation

    AttnInput = mean(Map(1).FX,2) + attnBaseInput; 
    attIdx = AttnInput > attnThreshold; 
    inhibVec(attIdx==1) = minInhib;  % rapid uncoupling when attention is directed to units
    inhibVec(attIdx==0) = P.spatinhib; % rapid re-coupling if attention is removed
    if (t > tlag*C.tstep), SpatInhib = sum(SpatAttnRecord(tIdx-tlag, inhibVec>0.5*P.spatinhib)); end   
    SpatAttn = min(maxAttn, max(0, SpatAttn + C.tstep * ( AttnInput - inhibVec.*SpatInhib) ));  % attraction of spatial attention to locations in feature maps, and global inhibition on spatial attention
    SpatAttnRecord(tIdx, :) = SpatAttn; 
    SpatInhibRecord(tIdx) = SpatInhib; 
    
%     subplot(1,2,1);
%     plot(AttnInput - inhibVec.*SpatInhib);
%     subplot(1,2,2);
%     plot(SpatAttn); 
%     
    % Map(1).FX = Map(1).FX + C.tstep * (-Map(1).FX + Map(1).FX .* repmat(SpatAttn, 1, C.nc)); % spatial attention modulates feature maps
    FStrength(tx, :) = sum(Map(1).FX(:, F(1:setsize)),1);
    SStrength(tx, :) = sum(SpatAttn);
    tx = tx + 1;
    
    if startNewC
        if currentCons > 0
            % ... encode the previous item now (one-shot encoding of the accumulated context and content)
            [B, GateClosed, gWeight, cIdx, newBindings] = IMencodeStim(B, context, content, GateClosed, gWeight, cRate(currentCons), cTime(currentCons));
            Strength(currentCons) = 0; % re-set, so that when this item is focused on again, it starts accumulating strength (in the FoA) over time from zero
            committed(currentCons, cIdx) = 1;
            bstrength(:, newBindings>0) = 0;  % set to 0 those Gating units that have been released (through delta) and re-recruited in the current encoding event
            bstrength(currentCons, :) = bstrength(currentCons, :) + newBindings;
            sumB = sum(bstrength(currentCons,:));
            lookAtMe = 1;
        end
        SAplus = SpatAttn + randn(C.nc, 1)*0.00000001; 
        PeakAttn = find(SAplus == max(SAplus));
        Distance2Items = abs(wrap(PeakAttn - L(1:setsize), 180));
        currentCons = find(Distance2Items == min(Distance2Items), 1);
        consSequence(cNum) = currentCons;
        if Inpos(currentCons) == 0, Inpos(currentCons) = inpos; inpos = inpos + 1; end   %keep record of the "input position", that is, the order in which items are (first) consolidated
        
        AfocusLoc = C.location(PeakAttn, :); % shift FoA to location with highest attentional weight in SpatAttn
        content = C.stimnoise;  % initialize new content
        if E.context == 1, context = C.locationnoise; end  % initialize new context
        if E.context == 2, context = C.stimnoise; end
        
        startNewC = 0;
        cClock = 0;
        cNum = cNum + 1;
    end
    
    if cClock >= cTime(cNum) % if the previous consolidation period is over, ...
        startNewC = 1;  % start new consolidation period
        SpatAttn = SpatAttn - P.IOR * C.location(PeakAttn, :)';
    end
    
    Afocus = AfocusLoc./sum(AfocusLoc) * Map(1).FX; % update what is read out from FX into the focus
    strength = (1-Strength(currentCons)) * (1-exp(-P.cRate*C.tstep)); % incremental strength: grows towards asymptote of 1
    Strength(currentCons) = Strength(currentCons) + strength;
    content = content + strength * Afocus * C.Mapping;
    if E.context == 1, context = context + strength * AfocusLoc * C.MappingC; end
    if E.context == 2, context = context + strength * (AfocusLoc./sum(AfocusLoc) * Map(2).FX) * C.Mapping; end
    
    tIdx = tIdx + 1;
    t = t + C.tstep;
    cClock = cClock + C.tstep;
end

if sum(Inpos==0) > 0
    Inpos(Inpos==0) = inpos:setsize;  % assign the items not consolidated at all the remaining input positions (in arbitrary order)
end

if plotTrajectory == 1
    disp([setsize, 0, consSequence]);
    
%     PreFigure([], [], 2);
%     for tIdx = 1:9
%        subplot(3,3,tIdx)
%        plot(SpatAttnRecord(5 + tIdx*10, :));
%     end
    
    PreFigure([], [], 2);
    subplot(1,2,1);
    unattendedL = find(attIdx==0); 
    unattendedL = unattendedL(randperm(length(unattendedL)));
    plot(C.tstep:C.tstep:E.RI, SpatAttnRecord(:, [unattendedL(1), L(1:setsize)]));
    subplot(1,2,2);
    plot(C.tstep:C.tstep:E.RI, SpatInhibRecord);
    
%     PreFigure([], [], 2);
%     subplot(1,2,1);
%     plot(C.tstep:C.tstep:E.RI, FStrength);
%     subplot(1,2,2);
%     plot(C.tstep:C.tstep:E.RI, SStrength);
end

Bstrength = mean(bstrength, 2);
finish = 1;






function [Map, L, F, CDAg, CDAw, Alpha, sumCtime, Pangle, Cangle, EEG_W, EEG_FX] = IMtrackSignals(setsize, Map, eW, eNoise)
% encodes a memory set simultaneously, simulates the proces time-step wise
% to track CDA and alpha over time

global P
global E
global C

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

P.asyFX = 5;
Timepoints = round(E.RI/C.tstep) + 1;
CDAg = zeros(1,Timepoints);
CDAw = zeros(1,Timepoints);
Alpha = zeros(1,Timepoints);
nElectrodes = size(eW, 2);
EEG_W = zeros(Timepoints, nElectrodes);
EEG_FX = zeros(Timepoints, nElectrodes);
context = rand(1, C.nCat)*0.1; % just to have a value at t=0

% determine whether a stimulus is integrated with a mask
maskWindow = ones(1, setsize) * P.maskWindow;
sdmaskWindow = P.maskWindowSD * P.maskWindow;
attentionWindows = gamrnd(maskWindow.^2./sdmaskWindow.^2, sdmaskWindow.^2./maskWindow, 1, setsize);
masking = attentionWindows > abs(E.MaskSOA);

strengthFX = randn(1, setsize) * P.SDstrengthFX + 1;

% Prepare Consolidation Times

cRate = gamrnd(P.cRate^2/P.cRateSD^2, P.cRateSD^2/P.cRate, 1, setsize+1);  % generate the vector of consolidation rates
cTime = -(log(1-P.cStrength)./cRate); % vector of consolidation times needed to reach strength P.cTau

if E.presentation == 1
    OnsetTimes = zeros(1, setsize+1);
    cTime = diff([0, min(E.RI, attentionWindows(1) + cumsum(cTime))]);
    TConsolidOnset = [attentionWindows(1), cumsum(cTime(1:setsize))];
    TConsolidEnd = [cumsum(cTime(1:setsize)), E.RI];
    FXupdated = [0, ones(1, setsize)]; % only for the first (and only) stimulus presentation an update of FX is triggered. When subsequent items are consolidated, no FX update occurs
end

if E.presentation == 2
    OnsetTimes = 0:(E.prestime + E.ISI):E.RI;
    OnsetTimes = OnsetTimes(1:setsize);
    ballistic = rand(1, setsize) < P.cBallistic;  % for each consolidation event determine whether it is ballistic or not

    % time lines for various events, time counting from the onset of the first stimulus
    TStim = [0:(E.prestime+E.ISI):((setsize-1)*(E.prestime+E.ISI)), setsize*(E.prestime+E.ISI) + E.RI];  % onset of stimuli; add the test event as last time because it curtails consolidation
    TattentionWindowEnd = TStim(1:setsize) + attentionWindows;
    TConsolidOnset = zeros(1, setsize);
    TConsolidEnd = zeros(1, setsize);
    TConsolidOnset(1) = TattentionWindowEnd(1);  % first consolidation ends after attention window of first stimulus
    TConsolidEnd(1) = TConsolidOnset(1) + cTime(1);
    if ballistic(1) == 0, TConsolidEnd(1) = min(TConsolidOnset(1) + cTime(1), E.prestime + E.ISI); end  % if not ballistic, cTime is cut short by next stimulus
    TestTime = setsize * (E.prestime + E.ISI) + E.RI;
    for inpos = 2:setsize
        if ballistic(inpos-1) == 1
            TConsolidOnset(inpos) = max(TStim(inpos), TConsolidEnd(inpos-1)); % ballistic consolidation: next consolidation can start only when previous has finished (and current stimulus has been shown)
            TConsolidEnd(inpos) = min(TConsolidOnset(inpos) + cTime(inpos), TestTime);
            if TConsolidOnset(inpos) > TStim(inpos + 1)
                TConsolidEnd(inpos) = TConsolidOnset(inpos);
            end
            % if by the time consolidation of item i can start, item i+1 has already been presented, the FX map will already have been updated,
            % and there is nothing to consolidate any more
        else
            TConsolidOnset(inpos) = TattentionWindowEnd(inpos);
            TConsolidEnd(inpos) = min(TConsolidOnset(inpos) + cTime(inpos), TStim(inpos+1));  % not ballistic: onset of next stimulus curtails consolidation
        end
    end
    TConsolidOnset = [TConsolidOnset, E.RI]; % take care of the case when inpos = setsize+1
    TConsolidEnd = [TConsolidEnd, E.RI]; % take care of the case when inpos = setsize+1
    cTime = max(0, TConsolidEnd - TConsolidOnset); % actually available consolidation time after accounting for consolidation postponenemt by preceding items, and curtailing by the test
    FXupdated = [zeros(1, setsize), 1];  % each time a new item is presented, FXupdated for that item is initially 0, so FX is updated (add a 1 at the end for when inpos > setsize)
end


sumCtime = sum(cTime);
consolStarted = zeros(1, setsize);
inpos = 1;
t = 0;
tcount = 1;


%%%%%%%%% Start Simulation Timestep by Timestep %%%%%%%%%%%%%%%%%%%555

while t < E.RI  % continue until end of RI
    % presentation of stimuli: encode into feature map
    if inpos <= setsize
        if t > OnsetTimes(inpos) && t <= (OnsetTimes(inpos) + E.prestime)  % in first iteration, t=0, so all measures pick up the pre-trial baseline
            if FXupdated(inpos) == 0
                Map = UpdateFX(Map);  %upon onset of the array, update FX (once!)
                FXupdated(inpos) = 1;
            end
            if E.presentation == 1
                % simultaneous presentation: parallel encoding into spatially organized feature maps;
                for ff = 1:C.nfeatures
                    maxFX = max(Map(ff).FX(:));
                    for item = 1:setsize
                        Map(ff).FX = Map(ff).FX + (P.asyFX - maxFX) * P.stimDrive * strengthFX(item) * C.location(L(item),:)' * (C.stim(F(ff,item),:));
                    end
                    % addition of mask if the mask falls within the replacement window
                    if masking(1) == 1  % masking is all-or-none for simultaneous array
                        for item = 1:setsize
                            Map(ff).FX = Map(ff).FX + (P.asyFX - maxFX) * P.stimDrive * C.location(L(item),:)' * C.maskStim;
                        end
                    end
                end
            end
            if E.presentation == 2
                maxFX = max(Map(ff).FX(:));
                item = inpos;
                Map(ff).FX = Map(ff).FX + (P.asyFX - maxFX) * P.stimDrive * strengthFX(item) * C.location(L(item),:)' * (C.stim(F(ff,item),:));
                % addition of mask if the mask falls within the replacement window
                if masking(item) == 1
                    Map(ff).FX = Map(ff).FX + (P.asyFX - maxFX) * P.stimDrive * C.location(L(item),:)' * C.maskStim;
                end
            end
        end
    end

    SpatAttn = mean(Map(1).FX,2);
    %SpatAttn = max( 0, SpatAttn + C.tstep * (mean(Map(1).FX,2) + P.TopDownSpatAttn.*AfocusLoc' - P.spatinhib*SpatAttn*sum(SpatAttn)) );  % attraction of spatial attention to locations in feature maps, top-down guidance by FoA, and global inhibition on spatial attention
    %Map(1).FX = Map(1).FX + C.tstep * (-Map(1).FX + Map(1).FX .* repmat(SpatAttn, 1, C.nc)); % spatial attention modulates feature maps

    if t >= TConsolidEnd(inpos)
        if inpos <= setsize
            % remove just-consolidated feature, so that anther one is the highest peak next (a form of inhibition of return)
            for ff = 1:E.nfeat
                Map(ff).FX = max(0, Map(ff).FX - P.IOR * AfocusLoc'*Afocus(ff,:));
            end
            if E.context == 2
                cueFromFX = (AfocusLoc./sum(AfocusLoc)) * Map(2).FX;
                Map(2).FX = max(0, Map(2).FX - P.IOR * AfocusLoc'*cueFromFX);
            end
        end
        inpos = inpos + 1; % move on to consolidation of next item
    end

    if t >= TConsolidOnset(inpos) && inpos <= setsize  % if the consolidation of item inpos starts, ...
        SpatAttn = mean(Map(1).FX,2);
        if consolStarted(inpos) == 0
            spatPeak = find(SpatAttn==max(SpatAttn), 1);    % find the peak of spatial attention ...
            AfocusLoc = C.ContextFun(C.x, deg2rad(spatPeak), P.kappaf_ctx);  % ... and move the FoA to that location
            sim2originalLoc = cosines(AfocusLoc', C.location(L(1:setsize), :)');
            Focus = find(sim2originalLoc == max(sim2originalLoc), 1);
            Afocus = AfocusLoc./sum(AfocusLoc) * Map(1).FX; % use location as (spatial) attentional filter to pull out the target feature from its feature map
            content = C.stimnoise + Afocus * C.Mapping;
            if E.context == 1, context = C.locationnoise + AfocusLoc * C.MappingC; end
            if E.context == 2, context = C.stimnoise + (AfocusLoc./sum(AfocusLoc) * Map(2).FX) * C.Mapping; end
            [W, G, GW, committedNew, ~] = IMencodeStim(W, context, content, G, GW, cRate(inpos), cTime(inpos), C.tstep);
            committedNewNotBase = randsample(committedNew, round(length(committedNew)*(1-P.pBase))); % sample a subset of newly committed binding units as the to-be-released ones
            consolStarted(inpos) = 1;  % set this to 1 so that the set-up of consolidation occurs only once per item
        end
    end
    if t > TConsolidOnset(1) % after consolidation of the first item has started, and some binding units have been committed ...
        % continuing release of binding units
        pLoss = 1 - 1/exp(P.rRate*C.tstep);  % see StepByStepEncoding.m
        nLoss = binornd(length(committedNewNotBase), pLoss);
        decommitted = randsample(committedNewNotBase, nLoss);
        committedNewNotBase = setdiff(committedNewNotBase, decommitted);
        G(decommitted) = 0;
        W(:, decommitted) = 0;  % remove weights to the now free binding units
    end

    [Map, W] = IMdecayFX(Map, W, C.tstep);   % decay of FX through one time step

    EEG_W(tcount,:) = (((context * W(1:C.nLocCat, :)) * W((C.nLocCat+1):end, :)') * C.Mapping') * eW + randn(1,nElectrodes)*eNoise;  % feed last-used context into weight matrix -> reactivate content -> project onto electrodes
    EEG_FX(tcount,:) = SpatAttn' * eW + randn(1,nElectrodes)*eNoise;
    CDAg(tcount) = sum(G);
    CDAw(tcount) = sum(abs(W(:)));
    Alpha(tcount) = sum(abs(EEG_FX(tcount,:)));
    tcount = tcount + 1;
    t = t + C.tstep;
end

Map = UpdateFX(Map);  % update Map during the test (which is not explicitly simulated here

Pangle = L(1:setsize);
Cangle = F(1:setsize);




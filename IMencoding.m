function [Map, W, GateClosed, gWeight, Focus, Afocus, content, context, Inpos, Strength, Bstrength, CTime, SpatAttn] = IMencoding(Map, W, GateClosed, gWeight, L, F, setsize, cueing, overTime)
% encodes a memory set, simultaneously or sequentially

if nargin < 9, overTime = 0; end %overTime is the time for ballistic consolidation that might have been taken by a preceding encoding event

global P
global E
global C

if nargin < 8, cueing = 0; end

% Pre-Cue:
CuedIdx = 0;
if (ismember(cueing, [2,3]) && E.PreRetro == 1)
    if (cueing == 2 || cueing == 5), CuedIdx = 1; end  % valid retro-cue: always on target
    if (cueing == 3)
        if (setsize == 1), CuedIdx = 1; else, CuedIdx = 1+randperm(setsize-1,1); end  % invalid retro-cue: never on target (except when setsize = 1)
    end
end

%  consolidation times
cRateSD = P.cRateSD*P.cRate;
cRate = gamrnd(P.cRate^2/cRateSD^2, cRateSD^2/P.cRate, 1, setsize);  % generate the vector of consolidation rates
cStrength = repmat(P.cStrength, 1, setsize);
if CuedIdx > 0
    cStrength(CuedIdx) = P.cStrength + (1-P.cStrength)*0.7; % pre-cued item receives higher cStrength -> longer consolidation
    cStrength(setdiff(1:setsize, CuedIdx)) = 0.5*P.cStrength; % not pre-cued items receive reduced cStrength
end

Inpos = zeros(1,setsize);  % initialize vector coding the input positions of the items (in the encoding order). In this vector, and the result vectors, the items are ordered by output position!
Focus = randperm(setsize, 1);  % start with random location of the focus
Afocus = zeros(E.nfeat,C.nc);            % just to initialize it
SpatAttn = zeros(C.nc, 1);
Strength = zeros(1,setsize);
committed = zeros(setsize, P.nb);   % keeps a record of which binding units are committed (i.e., gate closed) for each item
Bstrength = zeros(setsize, P.nb);

% determine whether a stimulus is integrated with a mask
maskWindow = ones(1, setsize) * P.maskWindow;
sdmaskWindow = P.maskWindowSD * P.maskWindow;
attentionWindows = gamrnd(maskWindow.^2./sdmaskWindow.^2, sdmaskWindow.^2./maskWindow, 1, setsize);
masking = attentionWindows > abs(E.MaskSOA);

strengthFX = randn(1, setsize) * P.SDstrengthFX + 1;

%%% Simultaneous presentation

if E.presentation == 1

    %if E.forwardrecall, encorder = 1:setsize; else, encorder = randperm(setsize); end

    Map = UpdateFX(Map);  % update Map in case this is a second array
    % simultaneous presentation: parallel encoding into spatially organized feature maps;
    % replacment by mask if the mask falls within the replacement window
    for ff = 1:C.nfeatures
        for inpos = 1:setsize
            Map(ff).FX = Map(ff).FX + strengthFX(inpos) * C.location(L(inpos),:)' * C.stim(F(ff,inpos),:);
        end
        P.asyFX = max(Map(ff).FX(:));
        if masking(1) == 1  % masking is all-or-none for simultaneous array
            %Map = UpdateFX(Map); % partial erasue of just-encoded stimuli
            for inpos = 1:setsize
                Map(ff).FX = Map(ff).FX + C.location(L(inpos),:)' * C.maskStim;
            end
        end
    end
    
    cTime = -(log(1-cStrength)./cRate); % vector of consolidation times needed to reach cStrength percent of maximal strength
    
    %cTime = diff([0, min(E.prestime + E.RI - max(attentionWindows(1), overTime), cumsum(cTime)) ]);
    % make sure that the cumulative consolidation time does not exceed presentation time + RI, minus the integration window because consolidation can start only after the integration window;
    % ballistic = rand < P.cBallistic;  % determine whether the consolidation event that hits E.RI is ballistic or not
    % if ballistic
    %     firstNotStarted = find([cTime, 0]==0, 1); % first consolidation event that does not even start because RI is exceeded
    %     extended = max(1, firstNotStarted-1); % --> give the preceding consolidation event (if there is one) an extension
    %     cTime(extended) = -log(1-P.cStrength)/cRate(extended);  % give that consolidation event the full time it needs to achieve the target strength
    % end

    ballistic = rand < P.cBallistic;  % determine whether consolidation of ALL items is ballistic or not
    if ballistic == 0
        cTime = diff([0, min(E.prestime + E.RI - max(attentionWindows(1), overTime), cumsum(cTime)) ]);
            % make sure that the cumulative consolidation time does not exceed presentation time + RI, minus the integration window because consolidation can start only after the integration window;
    end


    %sequential encoding - binding (read out from FX in sequential order)
    cumTime = 0;
    content = zeros(E.nfeat, C.nCat);
    context = zeros(1, C.nLocCat);
    outpos = randperm(setsize);   % random output order

    for inpos = 1:setsize

        if (cTime(inpos) > 0) && C.consolidAttempt >= inpos

            % focus on the highest spatial peak
            SpatAttn = mean(Map(1).FX,2);
            %SpatAttn = max(0, SpatAttn - P.spatinhib*SpatAttn*sum(SpatAttn));

            if max(SpatAttn) > 0

                % in case of pre-cue, start consolidating with this item
                if (ismember(cueing, [2,3]) && E.PreRetro == 1 && inpos == 1)
                    Focus = CuedIdx;
                    AfocusLoc = C.location(L(Focus),:);     % location attended to in the feature maps
                else
                    spatPeak = find(SpatAttn==max(SpatAttn), 1);
                    AfocusLoc = C.ContextFun(C.x, deg2rad(spatPeak), P.kappaf_ctx);
                    sim2originalLoc = cosines(AfocusLoc', C.location(L(1:setsize), :)');
                    Focus = find(sim2originalLoc == max(sim2originalLoc), 1);
                end
                %Focus = encorder(inpos);  % Focus codes the item index for the attended item, which is the (planned) output position (in case of full report)
                % focus attention to the relevant stimulus location in the feature map(s) FX
                %AfocusLoc = C.location(L(Focus),:);     % location attended to in the feature maps

                % retrieve contents from feature maps --> consolidation
                for ff = 1:E.nfeat
                    Afocus(ff,:) = AfocusLoc./sum(AfocusLoc) * Map(ff).FX;
                    content(ff,:) = C.stimnoise + Afocus(ff,:) * C.Mapping;
                end % use location as (spatial) attentional filter to pull out the target feature from its feature map

                % define the context (i.e., the anticipated retrieval cue) and the content
                if E.context == 1, context = C.locationnoise + AfocusLoc * C.MappingC; end  %
                if E.context == 2
                    cueFromFX = (AfocusLoc./sum(AfocusLoc)) * Map(2).FX;
                    context = C.stimnoise + cueFromFX * C.Mapping;
                end

                % bind content(s) and context together through mediated binding
                if C.seqVariant == 2 && inpos == setsize, rTime = max(cTime(inpos), E.prestime + E.RI - cumTime); else, rTime = cTime(inpos); end
                [W, GateClosed, gWeight, cIdx, bStrength] = IMencodeStim(W, context, content, GateClosed, gWeight, cRate(inpos), cTime(inpos), rTime); % ... now actually compute the result of consolidation in the preceding period
                committed(inpos, cIdx) = 1;
                Bstrength(1:(inpos-1), bStrength>0) = 0;
                Bstrength(inpos, :) = bStrength;
                Strength(outpos(inpos)) = 1-exp(-cRate(inpos).*cTime(inpos)); % effectively reached cStrength in light of the actually realized cTime

                % remove target feature, so that anther one is the highest peak next (a form of inhibition of return)
                for ff = 1:E.nfeat
                    Map(ff).FX = max(0, Map(ff).FX - P.IOR * AfocusLoc'*Afocus(ff,:));
                end
                if E.context == 2
                    cueFromFX = (AfocusLoc./sum(AfocusLoc)) * Map(2).FX;
                    Map(2).FX = max(0, Map(2).FX - P.IOR * AfocusLoc'*cueFromFX);
                end

                timePassing = cTime(inpos);
                cumTime = cumTime + timePassing;
                decaytime = min(max(0, cumTime-E.prestime), timePassing);
                [Map, W] = IMdecayFX(Map, W, decaytime);
                Inpos(outpos(inpos)) = inpos;  % the order of entries is the (anticipated) output order, which is random

            end

        end

    end

    Inpos(Inpos==0) = setdiff(1:setsize, Inpos); % fill empty Inpos values with the remaining, not consolidated items

    % pre-cued item is brought to the FoA after all items have been encoded
    if (ismember(cueing, [2,3]) && E.PreRetro == 1)
        Focus = CuedIdx;
        AfocusLoc = C.location(L(Focus),:);     % location attended to in the feature maps
        for ff = 1:E.nfeat, Afocus(ff,:) = AfocusLoc./sum(AfocusLoc) * Map(ff).FX; end % use location as (spatial) attentional filter to pull out the target feature from its feature map
    end

    decaytime = max(0, E.RI - max(0, cumTime-E.prestime)); % remaining time until test; cumTime-E.prestime is the part of the RI that has already been spent with consolidation
    [Map, W] = IMdecayFX(Map, W, decaytime);

end

%%% Sequential presentation

if E.presentation == 2
    if E.forwardrecall, encorder = 1:setsize; else, encorder = randperm(setsize); end
    ballistic = rand(1, setsize) < P.cBallistic;  % for each consolidation event determine whether it is ballistic or not

    if C.seqVariant == 1 || C.seqVariant == 2
        cTime = -(log(1-cStrength)./cRate); % vector of consolidation times needed to reach 90% of maximal strength

        % time lines for various events, time counting from the onset of the first stimulus
        TStim = [0:(E.prestime+E.ISI):((setsize-1)*(E.prestime+E.ISI)), setsize*(E.prestime+E.ISI) + E.RI];  % add the test event as last time because it curtails consolidation
        TattentionWindowEnd = TStim(1:setsize) + attentionWindows;
        TConsolidOnset = zeros(1, setsize);
        TConsolidEnd = zeros(1, setsize);

        TConsolidOnset(1) = TattentionWindowEnd(1);  % first consolidation ends after attention window of first stimulus
        TConsolidEnd(1) = TConsolidOnset(1) + cTime(1);
        if ballistic(1) == 0, TConsolidEnd(1) = min(TConsolidOnset(1) + cTime(1), E.prestime + E.ISI - overTime); end  % if not ballistic, cTime is cut short by next stimulus
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
        cTime = max(0, TConsolidEnd - TConsolidOnset); % actually available consolidation time after accounting for consolidation postponenemt by preceding items, and curtailing by the test
    end
    if C.seqVariant == 3
        cTime = max(0, repmat(E.prestime+E.ISI, 1, setsize) - attentionWindows);
    end
    if C.seqVariant == 1
        rTime = diff([TConsolidOnset, TestTime]); % if cTime exceeds prestime (i.e., consolidation bleeds into the presentation interval of the next stimulus presentation), so does release of BP units
    else
        rTime = max(0, repmat(E.prestime+E.ISI, 1, setsize) - attentionWindows);
    end

    cumTime = 0;

    for inpos = 1:setsize

        Focus = encorder(inpos); % to-be-encoded item (not necessarily = Focus!)
        if (inpos > 1), Map = UpdateFX(Map); end  % reset of feature map(s) because a new stimulus is attended
        for ff = 1:C.nfeatures
            Map(ff).FX = Map(ff).FX + strengthFX(inpos) * C.location(L(Focus),:)' * C.stim(F(ff,Focus),:);
            P.asyFX = max(Map(ff).FX(:));
            if masking(inpos) == 1  % masking is all-or-none for simultaneous array
                %Map = UpdateFX(Map); % partial erasue of just-encoded stimuli
                Map(ff).FX = Map(ff).FX + C.location(L(Focus),:)' * C.maskStim;
            end
        end

        SpatAttn = mean(Map(1).FX,2);
        %SpatAttn = max(0, SpatAttn - P.spatinhib*SpatAttn*sum(SpatAttn));

        % focus attention to the relevant stimulus location in the feature map(s) FX
        AfocusLoc = C.location(L(Focus),:);     % update location attended to in the feature maps

        % retrieve contents from feature maps --> consolidation
        for ff = 1:E.nfeat
            Afocus(ff,:) = AfocusLoc./sum(AfocusLoc) * Map(ff).FX;
            content(ff,:) = C.stimnoise + Afocus(ff,:) * C.Mapping;
        end % use location as (spatial) attentional filter to pull out the target feature from its feature map

        % define the context (i.e., the anticipated retrieval cue) and the content
        if E.context == 1, context = C.locationnoise + AfocusLoc * C.MappingC; end  %
        if E.context == 2
            cueFromFX = (AfocusLoc./sum(AfocusLoc)) * Map(2).FX;
            context = C.stimnoise + cueFromFX * C.Mapping;
        end

        if cTime(inpos) > 0
            [W, GateClosed, gWeight, cIdx, bStrength] = IMencodeStim(W, context, content, GateClosed, gWeight, cRate(inpos), cTime(inpos), rTime(inpos)); % ... now actually compute the result of consolidation in the preceding period
            committed(inpos, cIdx) = 1;
            Bstrength(1:(inpos-1), bStrength>0) = 0;
            Bstrength(inpos, :) = bStrength;
        end
        Inpos(Focus) = inpos;   %keep record of the input position
        Strength(Focus) = 1-exp(-cRate(inpos).*cTime(inpos)); % effectively reached cStrength in light of the actually realized cTime

        % remove target feature (a form of inhibition of return)
        for ff = 1:E.nfeat
            Map(ff).FX = max(0, Map(ff).FX - P.IOR * AfocusLoc'*Afocus(ff,:));
        end
        if E.context == 2
            cueFromFX = (AfocusLoc./sum(AfocusLoc)) * Map(2).FX;
            Map(2).FX = max(0, Map(2).FX - P.IOR * AfocusLoc'*cueFromFX);
        end

        decaytime = max(cTime(inpos), E.prestime+E.ISI-attentionWindows(inpos));
        if (inpos < setsize), [Map, W] = IMdecayFX(Map, W, decaytime); end % regardless of consolidation, decay time between one stimulus and the next entering FX is always ISI

    end

    [Map, W] = IMdecayFX(Map, W, E.RI);

end

KeepFocus = rand < P.keepFocus;
%if (setsize==1), KeepFocus = 1; end
Afocus = KeepFocus*Afocus;

CTime = mean(cTime);
Bstrength = mean(Bstrength, 2);


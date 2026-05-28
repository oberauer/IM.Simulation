function [Map, w, Focus, cTime, Inpos, Strength] = IMencoding(Map, w, L, F, setsize, Pgrad)
% encodes a memory set, simultaneously or sequentially

global P
global E
global C


cRate = gamrnd(P.cRate^2/P.cRateSD^2, P.cRateSD^2/P.cRate, 1,setsize);  % generate the vector of consolidation rates
cTime = -(log(1-P.cTau)./cRate); % vector of consolidation times needed to reach strength P.cTau
preMaskTime = min(cTime, E.MaskSOA); % time for consolidation before the mask

if E.presentation == 1  % simultaneous presentation: parallel encoding into spatially organized feature maps
    for ff = 1:C.nfeatures
        for inpos = 1:setsize
            Map(ff).FX = Map(ff).FX + C.location(L(inpos),:)' * C.stim(F(ff,inpos),:);
        end
    end
end

%sequential encoding - binding (read out from FX in sequential order)
Inpos = zeros(1,setsize);  % initialize vector coding the input positions of the items (in the encoding order).
% In this vector, and the result vectors, the items are ordered by output position!

Focus = randperm(setsize, 1);  % start with random location of the focus
AfocusLoc = C.location(L(Focus),:);     % location attended to
Afocus = zeros(1,C.nc);                 % just to initialize it
Strength = zeros(1,setsize);

if E.forwardrecall, encorder = 1:setsize; else, encorder = randperm(setsize); end % for consolidation: array items are encoded sequentially (in a random order chosen by the subject). This variable codes the planned output position of the encoded items

if E.presentation == 1
    for inpos = 1:setsize
        Focus = encorder(inpos);  % Focus codes the item index for the attended item, which is the (planned) output position (in case of full report)
        % focus attention to the relevant stimulus location in the feature map(s) FX
        AfocusLoc = C.location(L(Focus),:);     % location attended to in the feature maps
        Afocus = AfocusLoc./sum(AfocusLoc) * Map(1).FX; % use location as (spatial) attentional filter to pull out the target feature from its feature map
        % define the context (i.e., the anticipated retrieval cue) and the content
        Strength(inpos) = 1-exp(-cRate(inpos)*(min(preMaskTime(inpos), cTime(inpos))));
        if E.context == 1, context = Strength(inpos) * AfocusLoc * C.MappingC + C.locationnoise; end  % if 1 feature is bound to its spatial location
        if E.context == 2, context = Strength(inpos) * (AfocusLoc./sum(AfocusLoc) * Map(2).FX) * C.Mapping + C.stimnoise; end % use location as (spatial) attentional filter to pull out the context feature from its feature map
        content = Strength(inpos) * Afocus * C.Mapping + C.stimnoise;  % noise is added independent of strength because it arises in the conceptual layer into which FX is read out/mapped
        % update feature maps with mask
        if E.MaskSOA < sum(cTime(1:inpos))
            Map = UpdateFX(Map, AfocusLoc);
            for item = 1:setsize, Map(1).FX = Map(1).FX + C.location(L(item),:)' * C.maskStim; end % add mask to FX in all stimulus locations
        end
        for ff = 1:C.nfeatures, Map(ff).FX = exp(-P.decay*cTime(inpos))*Map(ff).FX; end  % decay
        % finally, bind content to context
        w = P.delta * w + (P.base + Pgrad*P.primacy.^(inpos-1)) .* ( context' * content );
        Inpos(Focus) = inpos;   %keep record of the input position
    end
end

if E.presentation == 2
    TStim = [0:E.prestime:((setsize-1)*E.prestime), inf];  % add inf at the end for the stimulus number setsize+1
    TMask = TStim + E.MaskSOA;
    cClock = 0; % clock for current consolidation process
    nextStim = 1;
    nextMask = 1;
    t = 0;
    inpos = 0;       % which item is currently presented
    currentCons = 0; % which item is currently being consolidated
    startNewC = 1;
    Strength = zeros(1,setsize);
    while t < min( max(sum(cTime), setsize*E.prestime), setsize*E.prestime + E.RI)  % continue until end of presentation, or end of all successive consolidations, whichever comes first
        % presentation of next stimulus: encode into feature map
        if t >= TStim(nextStim) % when the time for the next stimulus has arrived...
            inpos = inpos + 1; % move on to next input position
            encItem = encorder(inpos); % to-be-encoded item (not necessarily = Focus!)
            if (inpos > 1), Map = UpdateFX(Map, AfocusLoc); end  % complete reset of feature map because a new stimulus is attended (decrement parameter = 0)
            for ff = 1:C.nfeatures
                Map(ff).FX = Map(ff).FX + C.location(L(encItem),:)' * C.stim(F(ff,encItem),:);  % add new stimuli to feature map(s)
            end
            Afocus = AfocusLoc./sum(AfocusLoc) * Map(1).FX; % update what is read out from FX into the focus
            nextStim = min(setsize+1, nextStim + 1);
            Inpos(encItem) = inpos;   %keep record of the input position
        end
        % presentation of next mask: encode into feature map
        if t >= TMask(nextMask)
            Map = UpdateFX(Map, AfocusLoc);
            Map(1).FX = Map(1).FX + C.location(L(Focus),:)' * C.maskStim; % add mask to FX in the location of the just-presented stimulus
            nextMask = nextMask + 1;
            Afocus = AfocusLoc./sum(AfocusLoc) * Map(1).FX; % update what is read out from FX into the focus
        end
        if (inpos > 1)
            if currentCons < inpos && cClock >= cTime(currentCons)  % if we're not yet consolidating the last-presented item, and the previous consolidation period is over, ...
                startNewC = 1;  % start new consolidation period
            end
        end
        if startNewC
            currentCons = inpos; % set currently consolidated item to the last presented item
            cClock = 0;    % reset consolidation clock
            startNewC = 0; % switch off the need to start new consolidation period
            Focus = encorder(inpos);   % focus on new stimulus location
            % focus attention to the relevant stimulus location in the feature map(s) FX
            AfocusLoc = C.location(L(Focus),:);     % update location attended to in the feature maps
            Afocus = AfocusLoc./sum(AfocusLoc) * Map(1).FX; % use location as (spatial) attentional filter to pull out the target feature from its feature map
            content = C.stimnoise;
            if E.context == 1, context = C.locationnoise; end
            if E.context == 2, context = C.stimnoise; end
        end
        strength = (1-Strength(currentCons)) * (1-exp(-cRate(currentCons)*C.tstep)); % incremental strength: grows towards asymptote of 1
        Strength(currentCons) = Strength(currentCons) + strength;
        content = content + strength * Afocus * C.Mapping;
        if E.context == 1, context = context + strength * AfocusLoc * C.MappingC; end
        if E.context == 2, context = context + strength * (AfocusLoc./sum(AfocusLoc) * Map(2).FX) * C.Mapping; end
        %w = P.delta * w + (P.base + Pgrad*P.primacy.^(inpos-1)) .* ( context' * content );    % encoding of the final item
        w = (1/C.meanNsteps) * (-P.delta*w + (P.base + Pgrad*P.primacy.^(inpos-1)) .* ( context' * content )); 
        
        % better: content = C.stimnoise; content = content + strength*(); w = w + increment *(-(1-P.delta) + (Pbase... context' * content)) 
        
        t = t + C.tstep;
        cClock = cClock + C.tstep;
    end
    
end

end


function [Map, W, GateClosed, gWeight, Focus, Afocus, lastrefreshed, CuedIdx, Strength] = IMcueing(Map, W, GateClosed, gWeight, Strength, Focus, Afocus, L, F, setsize, cueing, refreshings)
% Does all the pre- and retro-cueing, refreshing, etc.

global P
global E
global C

lastrefreshed = 0;  %default
CuedIdx = 0; % default

if C.retroCueConsolid > 0
    % consolidation rate of the cued item
    cRate = gamrnd(P.cRate^2/P.cRateSD^2, P.cRateSD^2/P.cRate, 1, setsize);
    if C.retroCueConsolid == 1
        cTime = -(log(1-(P.cStrength-Strength))./cRate); % consolidation time needed to reach P.cStrength percent of maximal strength, taking into account the Strength already reached at encoding
        cTime = cTime .* (Strength==0);  % any item already consolidated, if only very briefly, is no longer eligible - otherwise, it gets full recruitment of additional binding units, with little release due to the usually very short cTime
    end
    if C.retroCueConsolid == 2
        cTime = -(log(1-(P.cStrength))./cRate); % consolidation time needed to reach P.cStrength percent of maximal strength, starting from zero because a second trace is laid down with a new set of binding units
    end
end

% determine whether the cue is integrated with FX or replaces FX
% integration is only possible for the first cue in a series - any
% subsequent cue must fall outside of the integration window

% sdIntegrate = P.integrateSD * P.integrate;
% integrationWindows = gamrnd(P.integrate.^2./sdIntegrate.^2, sdIntegrate.^2./P.integrate);
% integration = integrationWindows > E.RI;

if E.targetDim == 2  % if the target dimension is orientation, ...
    UpdateFX(Map, 1);  % the cue erases the target feature map, because that is the target feature map
end

if (ismember (cueing, [1:3, 5]))  % if the cueing condition is NOT "refreshing", and NOT multi-cueing

    %During retention interval - focus on a single item
    if (cueing == 1), CuedIdx = 0; end  % no cue - no cued item is yet defined
    if (cueing == 2 || cueing == 5), CuedIdx = 1; end  % valid retro-cue: always on target
    if (cueing == 3)
        if (setsize == 1), CuedIdx = 1; else, CuedIdx = 1+randperm(setsize-1,1); end  % invalid retro-cue: never on target (except when setsize = 1)
    end

    if CuedIdx > 0 && CuedIdx ~= Focus    % if a location was cued and forced the focus to shift during the RI...
        Focus = CuedIdx;          % ...update focus location, and ...
        Afocus = zeros(1,C.nc);   % ... drop content of the feature focus, as well as ...
        %AfocusLoc = C.location(L(Focus),:);  % ... update context (location) focus
    end

    if E.CTI(cueing) > 0
        [Afocus, AfocusLoc] = Retrieve(W, Map, Focus); 
        if C.retroCueConsolid > 0
            if cTime(Focus) > 0 % if the cued item has not yet been consolidated
                if E.context == 1, context = C.locationnoise + AfocusLoc * C.MappingC; end  %
                if E.context == 2, context = C.stimnoise + (AfocusLoc./sum(AfocusLoc) * Map(2).FX) * C.Mapping; end
                content = C.stimnoise + Afocus * C.Mapping;
                [W, GateClosed, gWeight, cIdx, bStrength] = IMencodeStim(W, context, content, GateClosed, gWeight, cRate(Focus), cTime(Focus), E.CTI(cueing));
            end
        end
        if P.cueingStrength > 0, W = Strengthen(W); end   % strengthening happens after retrieval, so it has no consequence for the cued item any more
        if P.removalThreshold > 0, [W, GateClosed] = Remove(W, GateClosed, setsize); end
    end

end

if (cueing == 4)    % "refreshing" experiments guiding the focus to several items during the retention interval

    %selection = randperm(length(C.RefSequence(refreshings).seq), 1);
    selection = randperm(size(C.RefSequence(refreshings).seq, 1), 1);
    refsequence = C.RefSequence(refreshings).seq(selection, :);  % select a refresh sequence
    for ref = 1:length(refsequence)
        Focus = refsequence(ref);   % focus goes to next item in the refresh sequence
        AfocusLoc = C.location(L(Focus),:); %  This entails an update of the location focus, and ...
        %Afocus = zeros(1,C.nc);              % ... the feature focus is dropped
        if C.retroCueConsolid > 0 && cTime(Focus) > 0  % if the cued item has not yet been consolidated
            for ff = 1:E.nfeat, Afocus(ff,:) = AfocusLoc./sum(AfocusLoc) * Map(ff).FX; end % use location as (spatial) attentional filter to pull out the target feature from its feature map
            if E.context == 1, context = C.locationnoise + AfocusLoc * C.MappingC; end  %
            if E.context == 2, context = C.stimnoise + (AfocusLoc./sum(AfocusLoc) * Map(2).FX) * C.Mapping; end
            for ff = 1:E.nfeat, content(ff,:) = C.stimnoise + Afocus(ff,:) * C.Mapping; end
            [W, GateClosed, gWeight, cIdx, bStrength] = IMencodeStim(W, context, content, GateClosed, gWeight, cRate(Focus), cTime(Focus), E.CTI(cueing));
            if C.retroCueConsolid == 1, cTime(Focus) = 0; end % now P.cStrength has certainly been reached, so no further consolidation is called for
        end
        if ref == length(refsequence), [Afocus, AfocusLoc] = Retrieve(W, Map, Focus); end  % no need to do that for pre-final refreshings because it has no consequence
        if P.cueingStrength > 0, W = Strengthen(W); end   % strengthening happens after retrieval, so it has no consequence for the last-cued item any more
        if P.removalThreshold > 0, [W, GateClosed] = Remove(W, GateClosed, setsize); end
        if ref < length(refsequence), [Map, W] = IMdecayFX(Map, W, E.CTI(cueing)); end
    end
    if (refreshings>1), lastrefreshed = find(refsequence==1, 1, 'last'); end
end

if (cueing == 6)    % multi-cueing (Rerko & Oberauer, 2013)

    CuedIdx = zeros(1,length(E.cuesequence));
    for cueNum = 1:length(E.cuesequence)
        if E.cuesequence(cueNum) == 0, CuedIdx(cueNum) = 1 + randperm(setsize-1, 1);  % select a non-target at random
        else, CuedIdx(cueNum) = E.cuesequence(cueNum); % select the target for the last cue
        end
        Focus = CuedIdx(cueNum);          % ...update focus location, and ...
        AfocusLoc = C.location(L(Focus),:); % ... update location focus as well, and
        %Afocus = zeros(1,C.nc);  % ... drop the feature focus
        if C.retroCueConsolid > 0 && cTime(Focus) > 0  % if the cued item has not yet been consolidated
            AfocusLoc = C.location(L(Focus),:);
            for ff = 1:E.nfeat, Afocus(ff,:) = AfocusLoc./sum(AfocusLoc) * Map(ff).FX; end % use location as (spatial) attentional filter to pull out the target feature from its feature map
            if E.context == 1, context = C.locationnoise + AfocusLoc * C.MappingC; end  %
            if E.context == 2, context = C.stimnoise + (AfocusLoc./sum(AfocusLoc) * Map(2).FX) * C.Mapping; end
            for ff = 1:E.nfeat, content(ff,:) = C.stimnoise + Afocus(ff,:) * C.Mapping; end
            [W, GateClosed, gWeight, cIdx, bStrength] = IMencodeStim(W, context, content, GateClosed, gWeight, cRate(Focus), cTime(Focus), E.CTI(cueing));
            if C.retroCueConsolid == 1, cTime(Focus) = 0; end % now P.cStrength has certainly been reached, so no further consolidation is called for
        end
        if cueNum == length(E.cuesequence), [Afocus, AfocusLoc] = Retrieve(W, Map, Focus); end  % no need to do that for pre-final refreshings because it has no consequence
        if P.cueingStrength > 0, W = Strengthen(W); end
        if P.removalThreshold > 0, [W, GateClosed] = Remove(W, GateClosed, setsize); end
        if cueNum < length(E.cuesequence), [Map, W] = IMdecayFX(Map, W, E.CTI(cueing)); end
    end

end


%%%%%%%%%%%% Embedded Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555

    function [Afocus, AfocusLoc] = Retrieve(wx, Map, Focus)
        % retrieve content at the cued context
        if E.context == 1
            AfocusLoc = C.location(L(Focus),:); % use currently focused location as retrieval cue (= CuedLoc)
            cue = [AfocusLoc * C.MappingC + C.locationnoise, zeros(1, C.nCat*E.nfeat)];
        end
        if E.context == 2   % for the case where the cue is a feature (rather than location)
            cue = [C.stim(F(2,Focus),:) * C.Mapping + C.stimnoise, zeros(1, C.nCat)]; % the 2nd feature is the retrieval cue for the first (= target) feature.
            AfocusLoc = ones(1,C.nc)*mean(C.location(L(Focus),:)); % uniformly distributed spatial focus with overall activation strength matching a locally focused focus
        end
        retrievedFX = AfocusLoc./sum(AfocusLoc) * Map(1).FX; % use location as (spatial) attentional filter for FX
        retrievedBinding = cue * wx;
        retrievedVec = retrievedBinding * wx';
        retrievedFeature = retrievedVec((C.nLocCat+1):(C.nLocCat+C.nCat));
        Afocus = retrievedFX + retrievedFeature * C.Mapping'; % vector of drift rates (one for each of the 360 colors) is computed as the strength with which each color is bound to the location cue

    end

    function wx = Strengthen(wx)

        % strengthening as retrieval and re-encoding
        % retrievedFX = AfocusLoc./sum(AfocusLoc) * Map(1).FX; % use location as (spatial) attentional filter for FX
        % Afocus = retrievedFX;
        % cue = [AfocusLoc * C.MappingC + C.locationnoise, zeros(1, C.nCat*E.nfeat)];
        % retrievedBinding = cue * wx;
        % retrievedVec = retrievedBinding * wx';
        % retrievedFeature = retrievedVec((C.nLocCat+1):(C.nLocCat+C.nCat));
        % Afocus = Afocus + retrievedFeature * C.Mapping';
        % content = C.stimnoise + Afocus * C.Mapping;
        % context = C.locationnoise + AfocusLoc * C.MappingC;
        % contentTransposed = content';
        % inputVec = [context, contentTransposed(:)'];   % the vectorization of the now vertical contents concatenates all features of multi-feature items
        % cuestrength = 1-exp(-P.cuerate*E.CTI(cueing));
        % wx = wx + P.cueingStrength .* cuestrength * (inputVec' * retrievedBinding);

        cuestrength = 1-exp(-P.cuerate*E.CTI(cueing));
        strengthenVec = [AfocusLoc * C.MappingC, zeros(1, E.nfeat*C.nCat)]; % the content side should be 0 so that nothing is added to wx
        strengthLoc = repmat(strengthenVec', 1, P.nb);  % strength with which each location category is strengthened
        wx = wx + P.cueingStrength*cuestrength * (strengthLoc .* wx);  %strengthening of bindings in focused location

    end

    function [wx, GateClosed] = Remove(wx, GateClosed, setsize)

        % removal
        % re-activate the binding units bound to the cued location
        cue = [AfocusLoc * C.MappingC + C.locationnoise, zeros(1, C.nCat*E.nfeat)];
        retrievedBinding = cue * W;
        removalStrength = logist((E.cuevalidity-1/setsize)./(1-1/setsize), P.removalTau, P.removalGain);
        releasedBinding = abs(retrievedBinding) < (removalStrength .* P.removalThreshold .* max(abs(retrievedBinding)));
        GateClosed(releasedBinding) = 0;
        wx(:, releasedBinding) = 0;
        % removing weights is included here for consistency with what happens in delta-updating. It reduces the strengthening benefit
        % and thereby makes it harder to find the sweet spot for the removal threshold (too high -> strengthening effect is destroyed;
        % too low --> removal benefit in reloading experiment disappears)
    end

end


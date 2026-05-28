function [response, rt, Map, B, Focus] = IMretrieve(Map, B, Focus, Afocus, probed, cueing, L, F, TestX, SameRange, ChangeRange, probeIdx)
% Retrieves one item

global E
global C
global P

A = zeros(1,C.nc);  % accumulator

% Use the last retro-cue, or the location last attended to at encoding, to build the retrieval cue
AfocusLoc = C.location(L(Focus),:); % use currently focused location as retrieval cue (= CuedLoc)
cue = AfocusLoc * C.MappingC + C.locationnoise;
Adrift = Afocus + (((cue * C.wcb) .* B) * C.wfb') * C.Mapping'; % vector of drift rates (one for each of the 360 colors) is computed as the strength with which each color is bound to the location cue

if E.test == 2  % if this is a recognition test
    YesNo = randn(1,2)*P.sz;  % initialize up the YesNo acccumulator for "same" vs. "change" decisions
end

t = 0;       % start the timer at 0 at the time point when the cue is presented
probePresented = 0;  % the probe has not yet been presented
proceed = 1;         % simulation should proceed

%%%% retrieval until the first feature reaches the boundary   %%%

while proceed   
    if t >= E.CTI(cueing) && probePresented == 0  % the cue-test interval has passed: now simulate what happens when the test display is presented
        probePresented = 1;  % register that the probe (and the color wheel) has been presented, so that the "if t > CTI" can be skipped next time
        Map = UpdateFX(Map);  % reset of feature map because a new stimulus (CW or probe) is attended
        Map(1).FX = Map(1).FX + TestX; %add test display with color wheel or probe to the target-feature map
        if Focus ~= probed   % if, at the onset of the probe, the focus is NOT on it, then ...
            Focus = probed;  % ... shift focus to the probe, and ...
            Afocus = zeros(1,C.nc); ... reset feature focus, and ...
            A = zeros(1,C.nc); % ... throw away all evidence accumulated so far, because it is misleading, and ...
        end
        if E.context == 1
            AfocusLoc = C.location(L(Focus),:); % set retrieval cue to the probed location
            cue = AfocusLoc * C.MappingC + C.locationnoise;
        end
        if E.context == 2
            cue = C.stim(F(2,Focus),:) * C.Mapping + C.stimnoise; % the 2nd feature is the retrieval cue for the first (= target) feature.
        end
        Afocus = Afocus + P.cwinter * AfocusLoc./sum(AfocusLoc) * Map(1).FX; % use location as (spatial) attentional filter for FX (FX contains only information from the test display)
        Afocus = Afocus + (((cue * C.wcb) .* B) * C.wfb') * C.Mapping';
        Adrift = Afocus ;     % compute a new drift rate
    end
    
    % Now we simulate the actual drift process: Each of the 360 evidence accumulators is incremented by a number of samples
    % that are drawn from a multinomial distribution.
    % The (normalized) drift rates are the probabilities with wich each accumulator receives an increment
    
    %A = A + (1+randn*P.driftnoise) * mnrnd(C.nsamples, Adrift./sum(Adrift))./C.nsamples; %sample from Adrift distribution, add with randomly varying drift rate
    %A = A + randn(1,C.nc) * P.dnoise; % add spontaneous noise to all accumulators
    
    A = A + (1+randn*P.driftnoise) * Adrift + randn(1,C.nc) * P.dnoise;  
    
    A = max(0, A - A*(P.inhib/C.nc)*sum(A)); % global inhibition (currently not used in the model because P.inhib = 0)
    
    t = t + C.tstep;   % increment time step for this simulation iteration
    %for ff = 1:C.nfeatures, Map(ff).FX = exp(-P.decay*C.tstep)*Map(ff).FX; end
    if max(A) > P.boundary(1) && probePresented, proceed = 0; end % check if boundary for retrieval has been reached
    if t > 5, proceed = 0; end  % if it takes too long, proceed after 5 s.
end

%%% If recall test: give the retrieved response %%%
if E.test == 1   % recall test
    response = find(A==max(A), 1);  % responses are ordered by output position ("probed" is incremented from 1 to E.outsize)
    rt = t-E.CTI(cueing);           % response time (subtract cue-target interval because t starts counting when the cue is presented, and RT should start when the target=probe is presented)
end

%%% If recognition test: Once question is given, accumulate Yes/No until 2nd boundary is reached %%%

if E.test == 2  %recognition: after deciding on which feature is retrieved, now decide on which response to give
    % after retrieval boundary has been reached: settle on one final retrieved feature value for comparison with the probe
    retrieved1 = find(A==max(A), 1);
    while max(YesNo) < P.boundary(2)  % while the decision boundary has not yet been reached
        if t > E.CQI(cueing)        % if the question has already been given
            v1 = [SameRange(retrieved1), ChangeRange(retrieved1)];   % now a unique feature value has actually been retrieved -> compare probe to that: Similar enough -> positive drift = 1 (otherwise 0); Dissimilar enough -> negative drift = 1 (otherwise 0)
            YesNo = YesNo + v1 + randn(1,2) * P.driftnoise;     % accumulate evidence towards the YES (=Same) and NO (=Change) responses
        end
        t = t + C.tstep;
    end  % exit the loop when response boundary has been reached
    response(1) = find(YesNo==max(YesNo), 1);  % now decide for response
    rt = t - max(E.CTI(cueing), E.CQI(cueing)); % record RT
    response(2) = abs(wrap(C.feature(F(1))-C.feature(probeIdx), 180)); % record the size of change of the probe relative to the target feature
    
end  % if recognition

B = P.delta * B; % output interference: reduce strength of binding pool activation

end


function [] = CreateMapping(calibrate)

global E
global C
global P

C.catCenter = linspace(0, C.nc, C.nCat+1) + randn(1, C.nCat+1)*P.mCatSD;  % feature at the centroid of each category
C.catCenter = max(1, min(C.nc, C.catCenter(1:C.nCat)));  % shave off the last element, which = C.nstim
KappaCateg = max(0.1, P.kappa_feat + randn(1,C.nCat)*P.kappaCatSD);  % create a set of kappas, one for each of the categories

% Mapping from sensory to categorical representation
% The mapping is also the distribution of features in each category!
C.Mapping = zeros(C.nc, C.nCat);
for cat = 1:C.nCat
    C.Mapping(:,cat) = VonMisesN(C.x, deg2rad(C.catCenter(cat)), KappaCateg(cat)); 
end

% create mapping from context stimuli to a more coarse (categorical) representation of context.
% No individual differences are assumed, so this can be done here.

LocCatCenter = linspace(0, C.nc, C.nLocCat+1);  % feature at the centroid of each category - for 8 categories (canonical orientations on the circle)
LocCatCenter = LocCatCenter(1:C.nLocCat);  % shave off the last element, which = C.nc

C.MappingC = zeros(C.nc, C.nLocCat);
for cat = 1:C.nLocCat
    C.MappingC(:,cat) = VonMisesN(C.x, deg2rad(LocCatCenter(cat)), P.kappa_ctx);
end

% Calibrate the mapping strength so that stimuli are reproduced with their
% original amplitude

if calibrate

    EContext = E.context;   % keep to restore later: Calibration always works with feature-location bindings (just to keep it simple)
    E.context = 1;

    nIterations = 500; % need to iterate many times to average out randomness in CreateMapping
    Amplify = zeros(1,nIterations);
    for iter = 1:nIterations
        C.amplify = 1;
        %CreateMapping;
        Tolerance = 0.01;
        Mismatch = 100;
        while (abs(Mismatch) > Tolerance)
            context = C.location(1,:) * (C.MappingC * C.amplify) + C.locationnoise;
            content = C.stim(round(C.nstim/2),:) * (C.Mapping * C.amplify) + C.stimnoise;

            % encoding
            W = CreateConnections(1);
            GateClosed = zeros(1, P.nb);  % gate-closing units
            GateWeight = zeros(1, P.nb); % gate-closing weights
            [W, ~, ~, ~, ~] = IMencodeStim(W, context, content, GateClosed, GateWeight, 1, 10, 10);

            % retrieval
            cue = [context, zeros(1, C.nCat)];
            retrievedBinding = cue * W;
            retrievedVec = retrievedBinding * W';
            retrievedFeature = retrievedVec((C.nLocCat+1):(C.nLocCat+C.nCat));

            retrieved = max(0, retrievedFeature * (C.Mapping * C.amplify)'); % vector of drift rates (one for each of the 360 colors) is computed as the strength with which each color is bound to the location cue

            Mismatch = (sum(C.stim(round(C.nstim/2),:)) - sum(retrieved))./sum(C.stim(round(C.nstim/2),:));  % proportional absolute deviation
            if Mismatch > 0
                C.amplify = C.amplify * 1.1;
            else
                C.amplify = C.amplify * 0.9;
            end
        end
        Amplify(iter) = C.amplify;
    end
    C.amplify = mean(Amplify);
    E.context = EContext;   %restore

end

C.Mapping = C.Mapping * C.amplify;
C.MappingC = C.MappingC * C.amplify;









function D = DualTaskConsolidation(featureOverlap)
% Simulation of Nieuwenstein & Wyble (2015)
% featureOverlap = 1: Decision task stimuli have feature overlap with
% memoranda (most experiments), = 0: no overlap (Experiment 5: auditory
% digits)

global P
global E
global C

E.mask = 1;
E.prestime = 0.1;
MaskSOA = [2, E.prestime];  %
setsize = 4;
E.outsize = setsize;
E.forwardrecall = 1; % participants were just asked to report all letters - presumably they did that mostly in forward order
C.consolidAttempt = setsize; % number of consolidation steps attempted in simultaneous condition
SOA = [0.25, 0.5, 1, 1.5];

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

Mdevobs = NaN(E.nsubj, 2, 2, length(SOA));  % id, singleDual, masking, III
MRT = NaN(E.nsubj, 2, length(SOA));  % id, masking, III
MCorrect = NaN(E.nsubj, 2, length(SOA));  % id, singleDual, masking, III
OT = struct('times', []); 
OverTime = repmat(OT, 1, length(SOA)); 

%[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x setsize x [1:360]

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;  % large number of colors to create the mask
    CreateMapping(E.calibrateAmp==2);
    while Cos > 0.2
        Stimulus = randn(2, C.nLocCat);  % stimuli for decision task
        Response = randn(2, C.nCat);     % responses for decision task
        CCS = cosines(Stimulus');
        CCR = cosines(Response');
        Cos = max(CCS(1,2), CCR(1,2));
    end

    for singleDual = 1:2
        for masking = 1:2
            for iii = 1:length(SOA)

                E.MaskSOA = MaskSOA(masking);
                if singleDual == 2
                    E.RI = SOA(iii) - E.prestime;
                else
                    E.RI = SOA(iii) - E.prestime + 1.2; % 1.2 s free time are added instead of the secondary task 
                end

                fdistance = zeros(E.ntrials, E.outsize);  % distance (target, response)
                rt = zeros(E.ntrials, 1);
                overTime = zeros(1, E.ntrials);
                correct = zeros(E.ntrials, 1);

                for trial = 1:E.ntrials

                    % map = struct('FX', zeros(C.nc));   % feature map
                    % Map = repmat(map, C.nfeatures, 1);
                    % W = CreateConnections(C.nfeatures);
                    % G = zeros(1, P.nb);  % gate-closing units
                    % GW = zeros(1, P.nb); % gate-closing weights
                    % 
                    % if singleDual == 2
                    %     % encode S-R bindings for decision task - doing this for every trial = loading the task set into WM
                    %     W2 = CreateConnections(C.nfeatures);
                    %     for stim = 1:2
                    %         [W2, G, GW] = IMencodeStim(W2, Stimulus(stim,:), Response(stim,:), G, GW, P.cRate, 1, 1); % consolidation time and release time = 1
                    %     end
                    % end
                    % 
                    % % memory array
                    % L(1,:) = 1:4;   % all in a row
                    % features = randperm(C.nstim);      %shuffle object features
                    % F = features(1:setsize);
                    % 
                    % [Map, W, G, GW, Focus, Afocus, content, context, Inpos, strength, bstrength, CTime, SpatAttn] = IMencoding(Map, W, G, GW, L, F, setsize, 1);
                    % usedTime = setsize*CTime; %CTime is the mean consolidation time taken
                    % overTime(trial) = max(0, usedTime - SOA(iii));
                    % 
                    % % decision task
                    % if singleDual == 2
                    %     selectedStim = randperm(2,1);
                    %     cue = [Stimulus(selectedStim,:), zeros(1, C.nCat)];
                    %     retrievedBinding = cue * W2;
                    %     retrievedVec = retrievedBinding * W2';
                    %     retrievedResponse = retrievedVec((C.nLocCat+1):(C.nLocCat+C.nCat));
                    %     Evidence = cosines(retrievedResponse', Response');
                    %     nsteps = round(5./C.tstep);   % that should be more than enough (5 sec)
                    %     drate = ones(nsteps, 1);
                    %     A = cumsum(drate * Evidence + randn(nsteps, 2)*P.dnoise);  % outer product of drate and Adrift -> matrix of tstep rows and 360 columns, each row = addition to to the 360 accumulators
                    %     maxA = max(A, [], 2); % maximum value in each row of A = maximum after each time step
                    %     t = find(maxA > P.boundary(2), 1);
                    %     if isempty(t), t=nsteps; end
                    %     response = find(A(t,:)==max(A(t,:)), 1);
                    %     rt(trial) = t.*C.tstep + overTime(trial);
                    %     correct(trial) = response == selectedStim;
                    % end
                    % 
                    % % memory test
                    % probed = 1;  % for now
                    % probestim = []; probeIdx = [];
                    % memresponse = IMretrieve(Map, W, G, Focus, Afocus, probed, 1, L, F, probestim, probeIdx);

                    output = IMSimDual(P, setsize, singleDual, Stimulus, Response, featureOverlap);  % featureOverlap = 1 indicates the feature overlap between memoranda and decision-task stimuli
                    for outpos = 1:E.outsize
                        fdistance(trial, outpos) = wrap(output.response(outpos)-output.F(outpos), 180);   %calculate distance between response and true feature in feature space (degrees!)
                    end
                    rt(trial) = output.drt;
                    correct(trial) = output.dcorrect;
                    overTime(trial) = output.overtime;

                end

                Mdevobs(id, singleDual, masking, iii) = mean(mean(abs(fdistance)));  %mean deviation (average over trials)
                if singleDual == 2
                    MRT(id, masking, iii) = mean(rt);
                    MCorrect(id, masking, iii) = mean(correct,1)';
                    OverTime(iii).times = [OverTime(iii).times, overTime];
                end

                disp('    ID        singDual  Masking   III       error     RT(dec.)   Acc(dec.)');
                disp([id, singleDual, masking, iii, mean(Mdevobs(id, singleDual, masking, iii)), (singleDual == 2).*MRT(id, masking, iii), (singleDual == 2).*MCorrect(id, masking, iii)]);

            end % III
        end  % masking
    end % single-dual

end  % for ID

MemAccuracy = (90-Mdevobs)./90; 

% Plot memory accuracy as function of III 
PreFigure;
for masking = 1:2
    subplot(1,2,masking);
    plotvector = squeeze(mean(MemAccuracy(:,:,masking,:),1));
    plot(SOA, plotvector);
    PostFigure([-0.1, max(SOA)+0.1, 0, 1], 'SOA', 'Accuracy', ['Masking = ', mat2str(masking-1)], {'Single', 'Dual'});
end

PreFigure
plotvector = squeeze(mean(MRT));
plot(SOA, plotvector);
PostFigure([-0.1, max(SOA)+0.1, 0, 1.2*max(plotvector(:))], 'SOA', 'RT(s)', [], {'No Mask', 'Mask'});

PreFigure;
for idx = 1:length(SOA)
    subplot(2,2,idx);
    h = histogram(OverTime(idx).times, 50);
    title('SOA = ', mat2str(SOA(idx)));
    text(0.5*h.BinLimits(2), 0.5*max(h.Values), ['Mean = ', mat2str(round(mean(OverTime(idx).times),5))]);
end

D.OverTime = OverTime;
D.MemAccuracy = MemAccuracy;
D.MRT = MRT;
D.MCorrect = MCorrect; 

halt = 1;
end

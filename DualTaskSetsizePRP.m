function D = DualTaskSetsizePRP(featureOverlap)
% Simulation of Experiment 1 of Stevanovsky & Jolicoeur (2007)
% featureOverlap = 1: Decision task stimuli have feature overlap with
% memoranda, = 0: no overlap 

global P
global E
global C

E.prestime = 0.1;
E.outsize = 1;
E.test = 2;   % change-detection
E.material = 2;
SOA = [0.1, 0.2, 0.4, 0.8];
Setsize = [1,2,4];
DecisionSetsize = [2,4];

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;
IMprepareRecog; % set up criterion for expected size of change

MemAccuracy = NaN(E.nsubj, 2, length(Setsize), length(DecisionSetsize), length(SOA));  % id, singleDual, setsize, decision-task setsize, III
MRT = NaN(E.nsubj, length(Setsize), length(DecisionSetsize), length(SOA));  % id, setsize, decision-task setsize, III
MCorrect = NaN(E.nsubj, length(Setsize), length(DecisionSetsize), length(SOA));  % id, setsize, decision-task setsize, III
MBindingResource = NaN(E.nsubj, length(Setsize), length(DecisionSetsize), length(SOA));  % id, setsize, decision-task setsize, III
OT = struct('times', []);
OverTime = repmat(OT, length(Setsize), length(SOA));

%[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x setsize x [1:360]

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;  % large number of colors to create the mask
    CreateMapping(E.calibrateAmp==2);

    for dssIdx = 1:length(DecisionSetsize)

        Cos = 1;
        while Cos > 0.2
            Stimulus = randn(DecisionSetsize(dssIdx), C.nLocCat);  % stimuli for decision task
            Response = randn(DecisionSetsize(dssIdx), C.nCat);     % responses for decision task
            CCS = cosines(Stimulus');
            CCR = cosines(Response');
            Cos = max(CCS(1,2), CCR(1,2));
        end

        for singleDual = 1:2
            for ssIdx = 1:length(Setsize)

                for soaIdx = 1:length(SOA)

                    setsize = Setsize(ssIdx);
                    if singleDual == 2
                        E.RI = SOA(soaIdx) - E.prestime;
                    else
                        E.RI = SOA(soaIdx) - E.prestime + 1.6; % 1.6 s free time are added instead of the secondary task
                    end

                    memcorrect = zeros(E.ntrials, 1);
                    rt = zeros(E.ntrials, 1);
                    overTime = zeros(1, E.ntrials);
                    correct = zeros(E.ntrials, 1);
                    G1 = zeros(E.ntrials, 1);

                    for trial = 1:E.ntrials

                        E.ptype = 1 + (rand > 0.5); % positive or new
                        output = IMSimDual(P, setsize, singleDual, Stimulus, Response, featureOverlap);
                        memresponse = output.response(1);
                        memcorrect(trial) = memresponse == E.ptype;
                        rt(trial) = output.drt;
                        correct(trial) = output.dcorrect;
                        overTime(trial) = output.overtime;
                        G1(trial) = mean(output.g1==0); % binding resource (from gating) free after encoding S-R mappings

                    end

                    MemAccuracy(id, singleDual, ssIdx, dssIdx, soaIdx) = mean(memcorrect);  %mean memory accoracy (average over trials)
                    if singleDual == 2
                        MRT(id, ssIdx, dssIdx, soaIdx) = mean(rt);
                        MCorrect(id, ssIdx, dssIdx, soaIdx) = mean(correct,1)';
                        OverTime(ssIdx, soaIdx).times = [OverTime(ssIdx, soaIdx).times, overTime];
                        MBindingResource(id, ssIdx, dssIdx, soaIdx) = mean(G1);
                    end

                    disp('    ID        singDual  setsize   D-setsize   III       error     RT(dec.)   Acc(dec.)');
                    disp([id, singleDual, Setsize(ssIdx), DecisionSetsize(dssIdx), SOA(soaIdx), mean(MemAccuracy(id, singleDual, ssIdx, dssIdx, soaIdx)), ...
                        (singleDual == 2).*MRT(id, ssIdx, dssIdx, soaIdx), (singleDual == 2).*MCorrect(id, ssIdx, dssIdx, soaIdx)]);

                end % III
            end  % setsize
        end % single-dual
    end % decision setsize

end  % for ID

% Plot memory accuracy as function of III

for dssIdx = 1:2
    PreFigure;
    for ssIdx = 1:length(Setsize)
        subplot(2,2,ssIdx);
        plotvector = squeeze(mean(MemAccuracy(:,:,ssIdx,dssIdx,:),1));
        plot(SOA, plotvector);
        PostFigure([0, max(SOA)+0.1, 0, 1], 'SOA', 'Accuracy', ['M-Setsize = ', mat2str(Setsize(ssIdx)), '; D-Setsize = ', mat2str(DecisionSetsize(dssIdx))], {'Single', 'Dual'});
    end
end

PreFigure
for dssIdx = 1:2
    subplot(1,2,dssIdx)
    plotvector = squeeze(mean(MRT(:,:,dssIdx,:)));
    plot(SOA, plotvector);
    PostFigure([0, max(SOA)+0.1, 0, 1.2*max(plotvector(:))], 'SOA', 'RT(s)', ['D-Setsize = ', mat2str(DecisionSetsize(dssIdx))], vec2legend(Setsize));
end

PreFigure
for dssIdx = 1:2
    subplot(1,2,dssIdx)
    plotvector = squeeze(mean(MCorrect(:,:,dssIdx,:)));
    plot(SOA, plotvector);
    PostFigure([0, max(SOA)+0.1, 0, 1.2*max(plotvector(:))], 'SOA', 'P(correct decision)', ['D-Setsize = ', mat2str(DecisionSetsize(dssIdx))], vec2legend(Setsize));
end

PreFigure
for dssIdx = 1:2
    subplot(1,2,dssIdx)
    plotvector = squeeze(mean(MBindingResource(:,:,dssIdx,:)));
    plot(SOA, plotvector);
    PostFigure([0, max(SOA)+0.1, 0, 1.2*max(plotvector(:))], 'SOA', 'Binging Units after Task Set', ['D-Setsize = ', mat2str(DecisionSetsize(dssIdx))], vec2legend(Setsize));
end

PreFigure;
for soaIdx = 1:length(SOA)
    OT = zeros(length(Setsize), length(OverTime(1,1).times));
    for ssIdx = 1:length(Setsize)
        OT(ssIdx,:) = OverTime(ssIdx,soaIdx).times;
    end
    subplot(2,2,soaIdx);
    hist(OT', 50);
    title('SOA = ', mat2str(SOA(soaIdx)));
end

Otime = zeros(length(Setsize), length(SOA));
for ssIdx = 1:length(Setsize)
    for soaIdx = 1:length(SOA)
        Otime(ssIdx, soaIdx) = mean(nonzeros(OverTime(ssIdx,soaIdx).times)); 
    end
end
PreFigure;
plot(SOA, Otime');
PostFigure([0, max(SOA)+0.1, 0, 1.2*max(Otime(:))], 'SOA', 'Overtime (s)', 'Ballistic Trials Only', vec2legend(Setsize));

D.OverTime = OverTime;
D.MemAccuracy = MemAccuracy;
D.MRT = MRT;
D.MCorrect = MCorrect; 
D.MBindingResource = MBindingResource;

halt = 1;
end

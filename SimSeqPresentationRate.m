function [] = SimSeqPresentationRate(Model, cRate, fitMM)
% Simulation of sequential vs. simultaneous presentation, continuous reproduction
% Varying the presentation rate for sequential presentation, and
% correspondingly, the RI for simultaneous presentation (as in Jacob's
% experiments)

global P
global E
global C

if ~exist('fitMM'), fitMM = 0; end

%E.mask = 1;
E.prestime = 0.15; 
E.MaskSOA = 0.15;  %
setsize = 4;
E.outsize = setsize;
P.cRate = cRate;
InterItemInterval = [0.2, 0.3, 0.4, 0.6, 1] - E.prestime;

option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

% initialize parameters of mixture model, IM, observed mean Deviation

MdevobsIn = NaN(E.nsubj, 2, length(InterItemInterval), setsize);  % id, sim-seq, III, inpos
MdevobsOut = NaN(E.nsubj, 2, length(InterItemInterval), setsize);  % id, sim-seq, III, outpos
Mdevobs1 = NaN(E.nsubj, 2, length(InterItemInterval), setsize);  % id, sim-seq, III, inpos
Mrt = NaN(E.nsubj, 2, length(InterItemInterval), setsize);  % id, sim-seq, III, inpos
MMSD = NaN(E.nsubj, 2, length(InterItemInterval));  % id, sim-seq, III
MMguessing = NaN(E.nsubj, 2, length(InterItemInterval));  % id, sim-seq, III
MMtranspos = NaN(E.nsubj, 2, length(InterItemInterval));  % id, sim-seq, III
Mwact = NaN(E.nsubj, 2, length(InterItemInterval));  % id, sim-seq, III
%SkippedC = NaN(E.nsubj, 2, length(InterItemInterval));  % id, sim-seq, III
meanCtime = NaN(E.nsubj, 2, length(InterItemInterval));  % id, sim-seq, III

[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x setsize x [1:360]

tcount = 1; %trial count

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;  % large number of colors to create the mask
    CreateMapping(E.calibrateAmp==2);

    for simseq = 1:2
        E.presentation = simseq;

        for iii = 1:length(InterItemInterval)   % 1 = simultaneous, 2+ = sequential

            if simseq == 2
                E.ISI = InterItemInterval(iii); E.RI = 0; 
            else
                E.RI = InterItemInterval(iii); 
                C.consolidAttempt = 4; % number of consolidation steps attempted in simultaneous condition
            end
            fdistance1 = NaN(E.ntrials, setsize);  % distance (target, response) for each input position (position in presentation order) in output position 1
            fdistanceIn = zeros(E.ntrials, setsize);  % distance (target, response) for each input position (position in presentation order) 
            fdistanceOut = zeros(E.ntrials, setsize);  % distance (target, response) for each output position (position in test order)
            RT = zeros(E.ntrials, setsize);       % RT for each input position
            wact = zeros(1,E.ntrials);
            Probedpos = zeros(E.ntrials,1);
            Pangle = zeros(E.ntrials,setsize);
            Cangle = zeros(E.ntrials,setsize+1);
            Targ = zeros(E.ntrials,1);
            Resp = zeros(E.ntrials,1);
            Setsize = zeros(E.ntrials,1);
            %SkippedConsolidation = zeros(E.ntrials,1);
            CTime = zeros(E.ntrials,1);

            for trial = 1:E.ntrials

                output = Model(P, setsize, 1);  % cueing = 1 (no cue)

                fdistance1(trial, output.Inpos(1)) =  wrap(output.response(1)-output.F(1,1), 180);  % output position 1 by input position
                for outpos = 1:setsize
                    fdistanceIn(trial, output.Inpos(outpos)) = wrap(output.response(outpos)-output.F(1,outpos), 180);   %calculate distance between response and true feature in feature space (degrees!)
                    fdistanceOut(trial, outpos) = wrap(output.response(outpos)-output.F(1,outpos), 180);   %calculate distance between response and true feature in feature space (degrees!)
                    RT(trial, output.Inpos(outpos)) = output.rt(outpos);
                    Strength(trial, output.Inpos(outpos)) = output.Strength(outpos); 
                end
                wact(trial) = sum(sum(output.g)); % sum of activation in weight matrix (in case of model 8, wfocus) -> CDA?
                tcount = tcount+1;

                %collect data for further modeling - only the first item tested
                Probedpos(trial) = output.L(1);
                Pangle(trial,:) = output.L(1:setsize);
                Cangle(trial,1:setsize) = output.F(1,1:setsize);
                Cangle(trial,setsize+1) = output.CWcolor;
                Targ(trial) = output.F(1, 1);
                Resp(trial) = output.response(1);
                Setsize(trial) = setsize;
                CTime(trial) = output.CTime; 
                %SkippedConsolidation(trial) = sum(output.Bstrength==0);
                %Bstrength(trial, :) = output.Bstrength';

            end


            MdevobsIn(id, simseq, iii, :) = mean(abs(fdistanceIn));  %mean deviation (average over output positions and trials)
            MdevobsOut(id, simseq, iii, :) = mean(abs(fdistanceOut));  %mean deviation (average over input position and trials)
            Mdevobs1(id, simseq, iii, :) = nanmean(abs(fdistance1));  %mean deviation in output position 1 (average over trials)

            Mrt(id, simseq, iii, :) = mean(RT);
            meanCtime(id, simseq, iii, :) = mean(CTime);
            %MBstrength(id, simseq, iii, :) = mean(Bstrength);

            Mwact(id, setsize) = mean(wact);

            ssData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting
            ssD.setsize = Setsize;
            ssD.response = Resp;
            ssD.Color = Cangle;

            if ~exist('Data', 'var')
                Data = ssData;
                D = ssD;
            else   % concatenate the data structures
                f = fieldnames(Data);
                for i = 1:length(f)
                    Data.(f{i}) = [Data.(f{i}); ssData.(f{i})];
                end
                ff = fieldnames(D);
                for i = 1:length(ff)
                    D.(ff{i}) = [D.(ff{i}); ssD.(ff{i})];
                end
            end

            % fit Mixture Model
            if fitMM
                startparms = [15, .1, .1, .1];
                lb = [eps, 0, 0, 0]; ub = [90, 1, 1, 1];
                npar = 3;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture, 4 -> Souza & Oberauer mixture (iC.ncluding color-wheel attraction)
                MMloglik = 500000;
                itercount = 0;
                while MMloglik > 400000
                    [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, ssData, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                    itercount = itercount + 1;
                end
                disp('    id        simseq    III       error     LL/1000   SD        P(unif)   P(swap)');
                disp([id, simseq, InterItemInterval(iii), mean(MdevobsIn(id, simseq, iii, :)), MMloglik/1000, MMparms]);
                MMSD(id, simseq, iii) = MMparms(1);
                MMguessing(id, simseq, iii) = MMparms(2);
                if npar > 2, MMtranspos(id, simseq, iii) = MMparms(3); end
            else
                disp('    id        simseq    constime  error');
                disp([id, simseq, iii, mean(MdevobsIn(id, simseq, iii, :))]);
            end

        end % III
    end % for simseq

end  % for ID

% Plot Mean(Deviation) as function of III 
% Stand-alone figures of errors, averaged over oputut positions 
PreFigure;
subplot(1,2,1);
plotvector = squeeze(mean(MdevobsIn(:,1,:,:),1));
meanvector = mean(plotvector, 2)'; 
plot(InterItemInterval, plotvector);
hold on
plot(InterItemInterval, meanvector, '-r'); 
PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 90], 'Inter-Item Interval', 'Deviation (Deg)', 'Sim. by Input Pos.', vec2legend(1:setsize));
subplot(1,2,2);
plotvector = squeeze(mean(MdevobsIn(:,2,:,:),1));
meanvector = mean(plotvector, 2)'; 
plot(InterItemInterval, plotvector);
hold on
plot(InterItemInterval, meanvector, '-r'); 
PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 90], 'Inter-Item Interval', 'Deviation (Deg)', 'Seq. by Input Pos', vec2legend(1:setsize));

% Plot Mean(Deviation) as function of III 
% Stand-alone figures of errors, averaged over input positions 
PreFigure;
subplot(1,2,1);
plotvector = squeeze(mean(MdevobsOut(:,1,:,:),1));
meanvector = mean(plotvector, 2)'; 
plot(InterItemInterval, plotvector);
hold on
plot(InterItemInterval, meanvector, '-r'); 
PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 90], 'Inter-Item Interval', 'Deviation (Deg)', 'Sim. by Output Pos', vec2legend(1:setsize));
subplot(1,2,2);
plotvector = squeeze(mean(MdevobsOut(:,2,:,:),1));
meanvector = mean(plotvector, 2)'; 
plot(InterItemInterval, plotvector);
hold on
plot(InterItemInterval, meanvector, '-r'); 
PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 90], 'Inter-Item Interval', 'Deviation (Deg)', 'Seq. By Output Pos', vec2legend(1:setsize));

% Plot Mean(Deviation) as function of III and input position - for output position 1 only
PreFigure;
subplot(1,2,1);
plotvector = squeeze(nanmean(Mdevobs1(:,1,:,:),1));
meanvector = mean(plotvector, 2)'; 
plot(InterItemInterval, plotvector);
hold on
plot(InterItemInterval, meanvector, '-r'); 
PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 90], 'Inter-Item Interval', 'Deviation (Deg)', 'Sim. by Input Pos. (Out=1)', vec2legend(1:setsize));
subplot(1,2,2);
plotvector = squeeze(nanmean(Mdevobs1(:,2,:,:),1));
meanvector = mean(plotvector, 2)'; 
plot(InterItemInterval, plotvector);
hold on
plot(InterItemInterval, meanvector, '-r'); 
PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 90], 'Inter-Item Interval', 'Deviation (Deg)', 'Seq. by Input Pos. (Out=1)', vec2legend(1:setsize));


PreFigure;
subplot(1,2,1);
plotvector = squeeze(mean(meanCtime(:,1,:),1));
plot(InterItemInterval, plotvector);
PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 1], 'Inter-Item Interval', 'Consol. Time/Item', 'Simultaneous');
subplot(1,2,2);
plotvector = squeeze(mean(meanCtime(:,2,:),1));
plot(InterItemInterval, plotvector);
PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 1], 'Inter-Item Interval', 'Consol. Time/Item', 'Sequential');


% Plot Mixture Model Parameters over Setsize
if fitMM
    MMPm = 1 - MMtranspos - MMguessing;
    %meanMMPm = squeeze(mean(MMPm,1));
    PreFigure;
    subplot(2,2,1);
    plot(InterItemInterval, squeeze(mean(MMSD(:,1,:),1)));
    PostFigure([-0.1, max(InterItemInterval)+0.1, 0, max(max(mean(MMSD,1)))+0.5], 'Inter-Item Interval', 'Mean SD', 'Simultaneous');
    subplot(2,2,2);
    plot(InterItemInterval, squeeze(mean(MMPm(:,1,:),1)));
    PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 1], 'Inter-Item Interval', 'Mean P(mem)', 'Simultaneous');
    subplot(2,2,3);
    plot(InterItemInterval, squeeze(mean(MMguessing(:,1,:),1)));
    PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 0.7], 'Inter-Item Interval', 'Mean P(guess)', 'Simultaneous');
    subplot(2,2,4);
    plot(InterItemInterval, squeeze(mean(MMtranspos(:,1,:),1)));
    PostFigure([0-0.1, max(InterItemInterval)+0.1, 0, 0.7], 'Inter-Item Interval', 'Mean P(trans)', 'Simultaneous');

    PreFigure;
    subplot(2,2,1);
    plot(InterItemInterval, squeeze(mean(MMSD(:,2,:),1)));
    PostFigure([-0.1, max(InterItemInterval)+0.1, 0, max(max(mean(MMSD,1)))+0.5], 'Inter-Item Interval', 'Mean SD', 'Sequential');
    subplot(2,2,2);
    plot(InterItemInterval, squeeze(mean(MMPm(:,2,:),1)));
    PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 1], 'Inter-Item Interval', 'Mean P(mem)', 'Sequential');
    subplot(2,2,3);
    plot(InterItemInterval, squeeze(mean(MMguessing(:,2,:),1)));
    PostFigure([-0.1, max(InterItemInterval)+0.1, 0, 0.7], 'Inter-Item Interval', 'Mean P(guess)', 'Sequential');
    subplot(2,2,4);
    plot(InterItemInterval, squeeze(mean(MMtranspos(:,2,:),1)));
    PostFigure([0-0.1, max(InterItemInterval)+0.1, 0, 0.7], 'Inter-Item Interval', 'Mean P(trans)', 'Sequential');

end


%%% Save results
if E.saveResults == 1
    fid = fopen('IMSim.SimSeqPresentationTime.dat', 'w');
    for id = 1:E.nsubj
        for simseq = 1:2
            for iii = 1:length(InterItemInterval)
                for inpos = 1:setsize
                    fprintf(fid, '%d %d %d  %d %d %d ', id, simseq, iii, inpos, MdevobsIn(id, simseq, iii, inpos), Mrt(id, simseq, iii, inpos));
                    %if (fitMM == 1), fprintf(fid, '%d %d %d', MMtranspos(id, simseq, iii), MMguessing(id, simseq, iii), MMSD(id, simseq, iii)); end
                    fprintf(fid, '\n');
                end
            end
        end
    end
    fclose(fid);
end


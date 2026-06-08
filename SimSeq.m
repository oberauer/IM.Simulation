function D = SimSeq(Model, fitMM)
% Simulation of sequential vs. simultaneous presentation, continuous reproduction
% Getting the effect right requires a rapid gate-closure rate (~20) so that
% sufficient information survives for consolidation after the mask in the
% simultaneous, separate condition.

global P
global E
global C

if ~exist('fitMM'), fitMM = 0; end

C.nfeatures = 2;
E.context = 2;  % use color to retrieve orientation, or vice versa
E.mask = 1;
E.MaskSOA = 0.08;  % 80 ms
E.prestime = E.MaskSOA;  
E.RI = 0.5 + 0.05; % including mask duration of 50 ms
setsize = 2;

option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

% initialize parameters of mixture model, IM, observed mean Deviation

Mdevobs = NaN(E.nsubj, 2, 2, 2);  % id, seq-sim, overlap, inpos
Mrt = NaN(E.nsubj, 2, 2);  % id, seq-sim, overlap
MMSD = NaN(E.nsubj, 2, 2);  % id, seq-sim, overlap
MMguessing = NaN(E.nsubj, 2, 2);  % id, seq-sim, overlap
MMtranspos = NaN(E.nsubj, 2, 2);  % id, seq-sim, overlap
Mwact = NaN(E.nsubj, 2, 2);  % id, seq-sim, overlap

ISICat = [0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8];  % lower bounds of ISI categories for sequential presentation
ISICatCount = zeros(E.nsubj, length(ISICat), 2);
MdevISI = zeros(E.nsubj, length(ISICat), 2); % id, ISI category, overlap

% For plotting response distributions:
% Array = zeros(E.nsubj*setsize*E.ntrials, setsize);
% Target = zeros(1,E.nsubj*setsize*E.ntrials);
% Response = zeros(1,E.nsubj*setsize*E.ntrials);

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

    for simseq = 1:2   % 1 = simultaneous, 2 = sequential
        for overlap = 1:2  % 1 = separate, 2 = overlapping in space

            E.presentation = simseq;
            E.layout = overlap;
            E.outsize = setsize;
            fdistance = zeros(E.ntrials, setsize, setsize);  % distance (target, response) for each input position (position in presentation order) x output position (position in test order)
            RT = zeros(E.ntrials, setsize, setsize);       % RT for each input x output position
            wact = zeros(1,E.ntrials);
            Probedpos = zeros(E.ntrials,1);
            Pangle = zeros(E.ntrials,setsize);
            Cangle = zeros(E.ntrials,setsize+1);
            Targ = zeros(E.ntrials,1);
            Resp = zeros(E.ntrials,1);
            Setsize = zeros(E.ntrials,1);
            featureDist = zeros(E.ntrials,2);

            for trial = 1:E.ntrials

                if simseq == 1 
                    isiCategory = 1;  % ISI = 0
                end
                if simseq == 2
                    E.ISI = 0.05 + rand * 0.65;  % 50 ms mask duration
                    isiCategory = find(E.ISI < ISICat, 1) - 1;
                end
                output = Model(P, setsize, 1);  % cueing = 1 (no cue)

                for outpos = 1:setsize
                    fdistance(trial, output.Inpos(outpos), outpos) = wrap(output.response(outpos)-output.F(1,outpos), 180);   %calculate distance between response and true feature in feature space (degrees!)
                    RT(trial, output.Inpos(outpos), outpos) = output.rt(outpos);
                end
                wact(trial) = sum(sum(output.g)); % sum of activation in weight matrix (in case of model 8, wfocus) -> CDA?
                tcount = tcount+1;

                %collect data for further modeling - only the first item tested
                Probedpos(trial) = output.F(2, 1);
                Pangle(trial,:) = output.F(2, 1:setsize);
                Cangle(trial,1:setsize) = output.F(1, 1:setsize);
                Cangle(trial,setsize+1) = output.CWcolor;
                Targ(trial) = output.F(1, 1);
                Resp(trial) = output.response(1);
                Setsize(trial) = setsize;

                featureDist(trial, :) = squeeze(sum(fdistance(trial,:,:), 2)); % sum over outpos, keep inpos
                MdevISI(id, isiCategory, overlap) = MdevISI(id, isiCategory, overlap) + mean(abs(featureDist(trial, :)));
                ISICatCount(id, isiCategory, overlap) = ISICatCount(id, isiCategory, overlap) + 1;

            end


            Mdevobs(id, simseq, overlap, :) = mean(abs(featureDist));  %mean deviation (averaged over trial)
            Mrt(id, simseq, overlap) = mean(nonzeros(RT));

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
                disp('    id        simseq    overlap   iteration error     LL/1000   SD        P(unif)   P(swap)');
                disp([id, simseq, overlap, itercount, mean(Mdevobs(id, simseq, overlap, :)), MMloglik/1000, MMparms]);
                MMSD(id, simseq, overlap) = MMparms(1);
                MMguessing(id, simseq, overlap) = MMparms(2);
                if npar > 2, MMtranspos(id, simseq, overlap) = MMparms(3); end
            else
                disp('    id        simseq    overlap   error');
                disp([id, simseq, overlap, mean(Mdevobs(id, simseq, overlap, :))]);
            end

        end  % for overlap
    end %for seqsim


end  % for ID

% Plot Mean(Deviation) as function of sequential-simultaneous, and of spatial overlap

% Stand-alone figures of errors
legendtext = {'Separate', 'Overlap'};
PreFigure;
subplot(1,2,1);
plotvector = squeeze(mean(Mdevobs(:,:,:,1),1));
plot(plotvector);
PostFigure([0.8, 2.2, 0, 1.05*max(max(plotvector))], 'Sim - Seq', 'Deviation (Deg)', 'First Presented', legendtext);
subplot(1,2,2);
plotvector = squeeze(mean(Mdevobs(:,:,:,2),1));
plot(plotvector);
PostFigure([0.8, 2.2, 0, 1.05*max(max(plotvector))], 'Sim - Seq', 'Deviation (Deg)', 'Second Presented', legendtext);

% plot errors for sequential presentation as a function of ISI category and
% overlap
sumDevISI = squeeze(sum(MdevISI,1));
sumISICatCount = squeeze(sum(ISICatCount,1));
plotvector = sumDevISI./sumISICatCount;
plotvector = plotvector(1:length(ISICat)-1,:); 
legendtext = {'Separate', 'Overlap'};
PreFigure;
plot(ISICat(1:length(ISICat)-1), plotvector);
PostFigure([-0.05, 1, 0, 1.05*max(max(plotvector))], 'ISI', 'Deviation (Deg)', [], legendtext);

D.Mdevobs = Mdevobs;
D.MDevISI = MDevISI;
D.ISICatCount = ISICatCount;

% Plot Mixture Model Parameters over Setsize
if fitMM
    MMPm = 1 - MMtranspos - MMguessing;
    %meanMMPm = squeeze(mean(MMPm,1));
    PreFigure;
    subplot(2,2,1);
    plot(squeeze(mean(MMSD,1)));
    PostFigure([0.8, 2.2, 0, max(max(mean(MMSD,1)))+0.5], 'Sim - Seq', 'Mean SD', 'SD from Bays Mixture', legendtext);
    subplot(2,2,2);
    plot(squeeze(mean(MMPm,1)));
    PostFigure([0.8, 2.2, 0, 1], 'Sim - Seq', 'Mean P(mem)', 'P(mem) from Bays Mixture', legendtext);
    subplot(2,2,3);
    plot(squeeze(mean(MMguessing,1)));
    PostFigure([0.8, 2.2, 0, 0.7], 'Sim - Seq', 'Mean P(guess)', 'P(guess) from Bays Mixture', legendtext);
    subplot(2,2,4);
    plot(squeeze(mean(MMtranspos,1)));
    PostFigure([0.8, 2.2, 0, 0.7], 'Sim - Seq', 'Mean P(trans)', 'P(trans) from Bays Mixture', legendtext);

    D.MMSD = MMSD;
    D.MMpm = MMPm;
    D.MMguessing = MMguessing;
    D.MMtranspos = MMtranspos;

end



%%% Save results
if E.saveResults == 1
    fid = fopen(['IMSim.SimSeq.dat'], 'w');
    for id = 1:E.nsubj
        for simseq = 1:2
            for overlap = 1:2
                for outpos = 1:setsize
                    fprintf(fid, '%d %d %d  %d %d  ', id, simseq, overlap, Mdevobs(id, simseq, overlap), Mrt(id, simseq, overlap));
                    if (fitMM == 1), fprintf(fid, '%d %d %d', MMtranspos(id, simseq, overlap), MMguessing(id, simseq, overlap), MMSD(id, simseq, overlap)); end
                    fprintf(fid, '\n');
                end
            end
        end
    end
    fclose(fid);
end


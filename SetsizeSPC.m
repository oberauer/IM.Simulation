function D = SetsizeSPC(Model, experiment, maxSetsize, fitMM)
% Simulation of Set-size and Serial-Position Effects for sequential
% presentation, continuous reproduction

global P
global E
global C

fitIM = 0;
fitIMSim = 0;

E.maxsetsize = maxSetsize;
E.presentation = 2;  % sequential!
E.prestime = 0.5*experiment; % Oberauer & Lin (2024): 0.5 s; Gorgoraptis et al. (2011): 0.5 + 0.5 s. 
E.context = experiment; % experiment 1: use location as cue; experiment 2: use other feature as cue
E.layout = experiment;  % 1: spatially distributed; 2: all in the same location
C.nfeatures = experiment;  % 1: only 1 feature; 2: need a 2nd feature as retrieval cue
% if (experiment == 1), P.keepFocus = 0.2; end  % for expectation of multiple tests
% if (experiment == 2), P.keepFocus = 0.8; end  % for expectation of a single test

if (experiment == 1), P.pBase = 0.5; end  % for expectation of multiple tests
if (experiment == 2), P.pBase = 0.2; end  % for expectation of a single test

option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

% initialize parameters of mixture model, IM, observed mean Deviation
Mdevobs = NaN(E.nsubj, E.maxsetsize, E.maxsetsize, E.maxsetsize);  % id, setsize, inpos, outpos
CircSD = NaN(E.nsubj, E.maxsetsize, E.maxsetsize, E.maxsetsize);  % id, setsize, inpos, outpos
Mrt = NaN(E.nsubj, E.maxsetsize, E.maxsetsize, E.maxsetsize);  % id, setsize, inpos, outpos
MMSD = NaN(E.nsubj, E.maxsetsize);  % id, setsize, inpos, outpos
MMguessing = NaN(E.nsubj, E.maxsetsize);  % id, setsize, inpos, outpos
MMtranspos = NaN(E.nsubj, E.maxsetsize);  % id, setsize, inpos, outpos
Mwact = NaN(E.nsubj, E.maxsetsize);  % id, setsize, inpos, outpos

IMparms = zeros(E.nsubj, 6);
IMSimparms = zeros(E.nsubj, 6);

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, E.maxsetsize), 1:360);  %Colors = E.ntrials x E.maxsetsize x [1:360]

tcount = 1; %trial count

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);

    for setsize = 1:E.maxsetsize

        if (experiment == 1), E.outsize = setsize; else, E.outsize = 1; end
        %fdistance = zeros(E.ntrials, setsize, E.outsize);  % distance (target, response) for each input position (position in presentation order) x output position (position in test order)
        fdistance = NaN(E.ntrials, setsize, E.outsize);  % distance (target, response) for each input position (position in presentation order) x output position (position in test order)
        RT = zeros(1,E.ntrials, setsize, E.outsize);       % RT for each input x output position
        wact = zeros(1,E.ntrials);
        Probedpos = zeros(E.ntrials,1);
        Pangle = zeros(E.ntrials,E.maxsetsize);
        Cangle = zeros(E.ntrials,E.maxsetsize+1);
        Targ = zeros(E.ntrials,1);
        Resp = zeros(E.ntrials,1);
        Setsize = zeros(E.ntrials,1);

        for trial = 1:E.ntrials

            output = Model(P, setsize, 1);  % cueing = 1 (no cue)

            for outpos = 1:E.outsize
                fdistance(trial, output.Inpos(outpos), outpos) = wrap(output.response(outpos)-output.F(outpos), 180);   %calculate distance between response and true feature in feature space (degrees!)
                RT(trial, output.Inpos(outpos), outpos) = output.rt(outpos);
            end
            wact(trial) = sum(sum(output.wx)); % sum of activation in weight matrix (in case of model 8, wfocus) -> CDA?
            tcount = tcount+1;

            %collect data for further modeling - only the first item tested
            if (E.context==1)
                Probedpos(trial) = output.L(1);
                Pangle(trial,:) = output.L(1:E.maxsetsize);
            end
            if (E.context==2)
                Probedpos(trial) = output.F(2,1);
                Pangle(trial,:) = output.F(2,1:E.maxsetsize);
            end
            Cangle(trial,1:setsize) = output.F(1,1:setsize);
            Cangle(trial,E.maxsetsize+1) = output.CWcolor;
            Targ(trial) = output.F(1);
            Resp(trial) = output.response(1);
            Setsize(trial) = setsize;

        end

        for inpos = 1:setsize
            for outpos = 1:E.outsize
                %Mdevobs(id, setsize, inpos, outpos) = mean(abs(fdistance(:,inpos,outpos)))*setsize;  %mean deviation - need to multiply by setsize because matrix is setsize*setsize but has only setsize non-zero entries
                Mdevobs(id, setsize, inpos, outpos) = nanmean(abs(fdistance(:,inpos,outpos)));  %mean deviation
                CircSD(id, setsize, inpos, outpos) = circ_std(deg2rad(fdistance(~isnan(fdistance(:,inpos,outpos)),inpos,outpos)));
                Mrt(id, setsize, inpos, outpos) = mean(RT(:,inpos,outpos))*setsize;
            end
        end
        Mwact(id, setsize) = mean(wact);

        ssData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting
        ssD.setsize = Setsize;
        ssD.response = Resp;
        if E.context == 1, ssD.L = round(C.Location(Pangle)); end
        if E.context == 2, ssD.L = round(Pangle); end
        ssD.Color = Cangle;

        if setsize == 1
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
            npar = 4;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture, 4 -> Souza & Oberauer mixture (iC.ncluding color-wheel attraction)
            MMloglik = 500000;
            itercount = 0;
            while MMloglik > 400000
                [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, ssData, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                itercount = itercount + 1;
            end
            disp([id, setsize, itercount, MMloglik]);
            MMSD(id, setsize) = MMparms(1);
            MMguessing(id, setsize) = MMparms(2);
            if npar > 2, MMtranspos(id, setsize) = MMparms(3); end
            if npar > 3, MMcwattraction(id, setsize) = MMparms(4); end
        end

        disp('     ID   setsize');
        disp([id, setsize]);

    end %for setsize

    %fit IM
    if fitIM
        startparms = [0.5, 0.5, 2, 10, 20, 0.5];  %B, A, s, P.kappa, P.kappafocus, Creduction
        npar = 6;
        lb = zeros(1,npar); ub = [5, 5, 20, 90, 90, 1];
        IMloglik = 500000;
        itercount = 0;
        while IMloglik > 400000
            [IMparms(id,:), IMloglik] = fminsearchbnd(@(x) IM(x, Data, 2), startparms, lb, ub, option);
            itercount = itercount + 1;
        end
        disp([id, itercount, IMloglik]);
        pred = IM(IMparms(id,:), Data, 1);
        Dev = abs(wrap(repmat(1:360, size(Data.response,1), 1) - repmat(Data.response, 1, 360), 180));
        predDev = sum(Dev .* pred, 2); % Weights deviation from each possible angle with that angle's predicted probability (of being the response). Computes weighted average deviation for predicted SD
        predMDevIM(id, :) = aggregate(Data.setsize, predDev);
    end

    %fit IMSim
    if fitIMSim
        startparms = [0.05, 1.5, 3, 10, 20, 0.5];  %X, Y, s, P.kappa, P.kappafocus, Creduction
        npar = 6;
        lb = zeros(1,npar); ub = [5, 5, 20, 90, 90, 1];
        IMSimloglik = 500000;
        itercount = 0;
        while IMSimloglik > 400000
            [IMSimparms(id,:), IMSimloglik] = fminsearchbnd(@(x) IMSim(x, D, 2), startparms, lb, ub, option);
            itercount = itercount + 1;
        end
        disp([id, itercount, IMSimloglik]);
        pred = IMSim(IMSimparms(id,:), D, 1);
        Dev = abs(wrap(repmat(1:360, size(D.response,1), 1) - repmat(D.Color(:,1), 1, 360), 180));
        predDev = sum(Dev .* pred, 2); % Weights deviation from each possible angle with that angle's predicted probability (of being the response). Computes weighted average deviation for predicted SD
        predMDevIMSim(id, :) = aggregate(D.setsize, predDev);

        pred2 = IMSim([X, Y, s, 2*P.kappa, 2*P.kappaf, r], D, 1);
        predDev2 = sum(Dev .* pred2, 2); % Weights deviation from each possible angle with that angle's predicted probability (of being the response). Computes weighted average deviation for predicted SD
        predMDevIMSim2(id, :) = aggregate(D.setsize, predDev2);

    end

end  % for ID

% Plot Mean(Deviation) as function of set size

% Stand-alone figures of errors, averaged over oputut positions (or input
% positions for plot of output position effect)
legendtext = {'In=1', 'In=2', 'In=3','In=4', 'In=5', 'In=6', 'In=7', 'In=8'};
PreFigure;
subplot(1,2,1);
plotvector = squeeze(nanmean(nanmean(Mdevobs,4),1));  % average over outpos (4) and subjects (1)
plot(plotvector);
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Deviation (Deg)', [], legendtext);

legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
subplot(1,2,2);
plot(plotvector');
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Input Position', 'Deviation (Deg)', [], legendtext);

if (experiment == 1)

    legendtext = {'Out=1', 'Out=2', 'Out=3','Out=4', 'Out=5', 'Out=6', 'Out=7', 'Out=8'};
    PreFigure;
    subplot(1,2,1);
    plotvector = squeeze(nanmean(nanmean(Mdevobs,3),1));  % average over inpos (3) and subjects (1)
    plot(plotvector);
    PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Deviation (Deg)', [], legendtext);

    legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
    subplot(1,2,2);
    plot(plotvector');
    PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Output Position', 'Deviation (Deg)', [], legendtext);

    PreFigure;
    for outpos = 1:setsize
        subplot(2,3,outpos);
        plotvector = squeeze(nanmean(Mdevobs(:,:,:,outpos),1))';
        plot(plotvector);
        PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Input Position', 'Deviation (Deg)', ['Output Position = ', mat2str(outpos)]);
    end

end

if (experiment == 2)

    PreFigure;
    subplot(1,2,1);
    plotvector = squeeze(nanmean(nanmean(Mdevobs(:,:,:,1),3),1));  % average over inpos (3) and subjects (1)
    plot(plotvector);
    PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Set Size', 'Deviation (Deg)');

    legendtext = {'In=1', 'In=2', 'In=3','In=4', 'In=5', 'In=6', 'In=7', 'In=8'};
    subplot(1,2,2);
    plotvector = squeeze(nanmean(nanmean(rad2deg(CircSD),4),1));  % average over outpos (4) and subjects (1)
    plot(plotvector);
    PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Circular SD (Deg)', [], legendtext);

    %     legendtext = {'In=1', 'In=2', 'In=3','In=4', 'In=5', 'In=6', 'In=7', 'In=8'};
    %     PreFigure;
    %     plotvector = squeeze(nanmean(nanmean(1./(CircSD),4),1));  % average over outpos (4) and subjects (1)
    %     plot(plotvector);
    %     PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Precision (1/rad)', [], legendtext);

    precX = NaN(E.maxsetsize, E.maxsetsize);
    for ss = 1:E.maxsetsize
        precX((E.maxsetsize-ss+1):E.maxsetsize, ss) = squeeze(mean(1./CircSD(:,ss,1:ss,1), 1));   % average over subjects (1)
    end
    legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
    PreFigure;
    plot(-E.maxsetsize:-1, precX);
    PostFigure([-(0.5+E.maxsetsize),-0.5, 0, 1.05*max(max(precX))], 'Recency', 'Precision(1/rad)', [], legendtext);
    % From circSD to precision is a non-linear transformation, changing the
    % order of set sizes with regard to the last-presented item: For larger
    % set sizes, circSD has more variability, leading to larger
    % mean(circSD) but also larger precision.

end

% Plot Mixture Model Parameters over Setsize
if fitMM
    MMPm = 1 - MMtranspos - MMguessing - MMcwattraction;
    meanMMPm = squeeze(mean(MMPm,1));
    meanK = bsxfun(@times, meanMMPm, 1:setsize);
    PreFigure
    subplot(3,2,1);
    plot(squeeze(mean(MMSD,1))');
    PostFigure([0.8,setsize+0.2, 0, max(max(mean(MMSD,1)))+0.5], 'Setsize', 'Mean SD', 'SD from Bays Mixture', {'Neutral', 'Valid', 'Invalid'});
    subplot(3,2,2);
    plot(meanK');
    PostFigure([0.8,setsize+0.2, 0, 6], 'Setsize', 'K', 'K from Bays Mixture', {'Neutral', 'Valid', 'Invalid'});
    subplot(3,2,3);
    plot(squeeze(mean(MMPm,1))');
    PostFigure([0.8,setsize+0.2, 0, 1], 'Setsize', 'Mean P(m)', 'P(mem)');
    subplot(3,2,4);
    plot(squeeze(mean(MMguessing,1))');
    PostFigure([0.8,setsize+0.2, 0, 0.5], 'Setsize', 'Mean P(guess)', 'P(guess)');
    subplot(3,2,5);
    plot(squeeze(mean(MMtranspos,1))');
    PostFigure([0.8,setsize+0.2, 0, 1], 'Setsize', 'Mean P(trans)', 'P(transpos)');
    subplot(3,2,6);
    plot(squeeze(mean(MMcwattraction,1))');
    PostFigure([0.8,setsize+0.2, 0, 1], 'Setsize', 'Mean P(wheel)', 'P(wheel attraction)');

    D.MMSD = MMSD;
    D.MMpm = MMPm;
    D.MMguessing = MMguessing;
    D.MMtranspos = MMtranspos;
    D.MMcwattraction = MMcwattraction;
end

% display parameter estimates in command window
if fitIM
    disp('      b       a          s       P.kappa      P.kappaf      r');
    disp(mean(IMparms, 1));
    disp('      b       a          s       P.kappa      P.kappaf      r');
    disp(std(IMparms, 1));
end

if fitIMSim
    disp('      X       Y          s       P.kappa      P.kappaf      r');
    disp(mean(IMSimparms, 1));
    disp('      X       Y          s       P.kappa      P.kappaf      r');
    disp(std(IMSimparms, 1));
end


%%% Save results
if E.saveResults == 1
    fid = fopen(['IMSim.SetsizeSPC', mat2str(experiment), '.dat'], 'w');
    for id = 1:E.nsubj
        for setsize = 1:E.maxsetsize
            if experiment==1, outsize = setsize; else, outsize = 1; end
            for inpos = 1:setsize
                for outpos = 1:outsize
                    fprintf(fid, '%d %d %d %d  %d %d %d  ', id, setsize, inpos, outpos, Mdevobs(id, setsize, inpos, outpos), CircSD(id, setsize, inpos, outpos), Mrt(id, setsize, inpos, outpos));
                    if (fitMM == 1), fprintf(fid, '%d %d %d %d', MMtranspos(id, setsize, inpos, outpos), MMguessing(id, setsize, inpos, outpos), MMcwattraction(id, setsize, inpos, outpos), MMSD(id, setsize, inpos, outpos)); end
                    fprintf(fid, '\n');
                end
            end
        end
    end
    fclose(fid);
end



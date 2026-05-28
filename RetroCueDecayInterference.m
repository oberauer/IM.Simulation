function [] = RetroCueDecayInterference(Model, setsize, fitMM, fitIMSim)
% Simulates the experiment of Hautekiet & Oberauer (2026) varying the RI of
% the no-cue condition, and eliminating perceptual interference (here:
% wheel = 0)

global C
global E
global P

option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);
E.cuevalidity = 1;
E.wheel = 0; % no percpetual interference at test
C.nloc = 5;  % only 5 equally spaced locations were used
RI = [1, 3, 5]; % levels of RI for no-cue condition

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

% initialize parameters of mixture model, IM, observed mean Deviation

Mdevobs =  NaN(E.nsubj, 2, 3); % Mean deviation (target feature, response) for subjects, cueing condition (no/valid), and RI
MMSD = NaN(E.nsubj, 2, 3);    % SD parameter of mixture model
MMguessing = NaN(E.nsubj, 2, 3);  % P(guess) parameter of mixture model
MMtranspos = NaN(E.nsubj, 2, 3);  % P(transposition) parameter of mixture model
MMcwattraction = NaN(E.nsubj, 2, 3);  % P(color wheel attraction) parameter of mixture model
MRT = NaN(E.nsubj, 2, 3);   % subj, cueing condition, wheel condition

IMSimparms = zeros(E.nsubj, 6);

Array = NaN(E.nsubj*2*3*E.ntrials, setsize+1);
Target = NaN(1,E.nsubj*2*3*E.ntrials);
Response = NaN(1,E.nsubj*2*3*E.ntrials);


[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x maxsetsize x [1:360]

tcount = 1; %trial count

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);

    for ri = 1:length(RI)   % 1: color wheel, 2: grey wheel

        fdistance = zeros(1,E.ntrials);
        rt = zeros(1,E.ntrials);
        Probedpos = zeros(E.ntrials,1);
        Pangle = zeros(E.ntrials,setsize);
        Cangle = zeros(E.ntrials,setsize);
        Targ = zeros(E.ntrials,1);
        Resp = zeros(E.ntrials,1);
        Setsize = zeros(E.ntrials,1);

        for cueing = 1:2  % 1: no cue (0 s delay), 2: valid cue (1 s delay),

            if cueing == 1, E.RI = RI(ri); end
            if cueing == 2, E.RI = 1; E.CSI(2) = RI(ri) - 1; end  % the CSI includes the 1 s presentation time of the cue and the post-cue time
            if cueing == 2 && ri == 1, ntrials = 0; else, ntrials = E.ntrials; end % a retro-cue condition with CSI=0 did not exist

            for trial = 1:ntrials

                output = Model(P, setsize, cueing);

                Array(tcount,:) = [output.F(1,1:setsize), output.CWcolor];  %add the color in the color wheel closest to target location in column setsize+1
                fdistance(trial) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
                rt(trial) = output.rt;
                Target(tcount) = output.F(1,1);
                Response(tcount) = output.response;
                tcount = tcount+1;

                %collect data for further modeling
                Probedpos(trial) = output.L(1);
                Pangle(trial,:) = output.L(1:setsize);
                Cangle(trial,1:setsize) = output.F(1,1:setsize);
                Cangle(trial,setsize+1) = output.CWcolor;
                Targ(trial) = output.F(1,1);
                Resp(trial) = output.response;
                Setsize(trial) = setsize;

            end

            if ntrials > 0

                Mdevobs(id, cueing, ri) = mean(abs(fdistance));  %mean deviation
                MRT(id, cueing, ri) = mean(rt);

                condData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for mixture-model fitting (separate structures for each condition)
                condData.cueing = repmat(cueing, length(Response), 1);
                condData.preretro = 2;
                condData.ri = ri;
                condD.setsize = Setsize;
                condD.response = Resp;
                condD.L = Pangle;
                condD.Color = Cangle;
                condD.cueing = repmat(cueing, length(Response), 1);
                condD.preretro = 2;
                condD.ri = ri;

                % Concatenate the condition-specific data to a data structure
                % across all conditions (for IM modelling)
                if ri == 1 && cueing == 1
                    D = condD;
                else   % concatenate the data structures
                    ff = fieldnames(D);
                    for i = 1:length(ff)
                        D.(ff{i}) = [D.(ff{i}); condD.(ff{i})];
                    end
                end

                % fit Mixture Model
                if fitMM
                    startparms = [15, .1, .1, .1];
                    lb = [eps, 0, 0, 0]; ub = [90, 1, 1, 1];
                    npar = 4;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture, 4 -> Souza & Oberauer mixture (including colorwheel attraction)
                    MMloglik = 500000;
                    itercount = 0;
                    while MMloglik > 400000
                        [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, condData, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                        itercount = itercount + 1;
                    end
                    disp('    ID        Setsize    Cueing   iter      LL/1000   SD        Pmem      Ptrans    Pwheel ');
                    disp([id, setsize, cueing, itercount, MMloglik/1000, MMparms(1), 1-sum(MMparms(2:4)), MMparms(3), MMparms(4)]);
                    MMSD(id, cueing, ri) = MMparms(1);
                    MMguessing(id, cueing, ri) = MMparms(2);
                    if npar > 2, MMtranspos(id, cueing, ri) = MMparms(3); end
                    if npar > 3, MMcwattraction(id, cueing, ri) = MMparms(4); end
                end

                disp('    id  cueing  RI ');
                disp([id, cueing, RI(ri)]);

            end

        end % for cueing
    end %for ri


    %fit IMSim
    if fitIMSim
        startparms = [0.05, 1.5, 3, 10, 20, 0.5];  %X, Y, s, kappa, kappafocus, Creduction
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

        pred2 = IMSim([X, Y, s, 2*kappa, 2*kappaf, r], D, 1);
        predDev2 = sum(Dev .* pred2, 2); % Weights deviation from each possible angle with that angle's predicted probability (of being the response). Computes weighted average deviation for predicted SD
        predMDevIMSim2(id, :) = aggregate(D.setsize, predDev2);

    end

end  % for ID

% Plot Mean(Deviation) as functions of Wheel
legendtext = {'No Cue', 'Retro-Cue'};

PreFigure
plotvector = squeeze(mean(Mdevobs,1));
plot(RI, plotvector);
PostFigure([0.8, max(RI)+0.2, 0, 1.05*max(max(plotvector))], 'RI', 'Deviation', 'Mean Deviation', legendtext);
% 
% 
% % Plot Mixture Model Parameters over Wheel condition
% if fitMM
%     MMPm = 1 - MMtranspos - MMguessing - MMcwattraction;
%     Xticks = {'Color', ' ', 'Grey'};  % Matlab insists on placing a tick every 0.5 steps, so have to insert a blank tick
%     PreFigure
%     subplot(3,2,1);
%     plot(squeeze(mean(MMSD,1))');
%     set(gca,'XtickLabel', Xticks);
%     PostFigure([0.8,2.2, 0, max(max(mean(MMSD,1)))+0.5], 'Wheel Condition', 'Mean SD', 'SD from Bays Mixture', {'No Cue', 'Retro-Cue'});
%     subplot(3,2,2);
%     plot(squeeze(mean(MMPm,1))');
%     set(gca,'XtickLabel', Xticks);
%     PostFigure([0.8,2.2, 0, 1], 'Wheel Condition', 'Mean P(m)', '"P(mem) Bays Mixture');
%     subplot(3,2,3);
%     plot(squeeze(mean(MMguessing,1))');
%     set(gca,'XtickLabel', Xticks);
%     PostFigure([0.8,2.2, 0, 0.5], 'Wheel Condition', 'Mean P(guess)', 'P(guess) Bays Mixture');
%     subplot(3,2,4);
%     plot(squeeze(mean(MMtranspos,1))');
%     set(gca,'XtickLabel', Xticks);
%     PostFigure([0.8,2.2, 0, 0.5], 'Wheel Condition', 'Mean P(trans)', 'P(transpos) Bays Mixture');
%     subplot(3,2,5);
%     plot(squeeze(mean(MMcwattraction,1))');
%     PostFigure([0.8,2.2, 0, 0.5], 'Wheel Condition', 'Mean P(CW)', 'P(CW) Bays Mixture');
% 
%     % Figure layout as in Souza et al (2016)
%     meanPm = squeeze(mean(MMPm, 1));
%     meanCWA = squeeze(mean(MMcwattraction, 1));
%     Xticks = {'No Cue', ' ', 'Cue'};
%     PreFigure;
%     subplot(2,2,1);
%     plot(1:2, meanPm(:,1)); % color-wheel
%     set(gca,'XtickLabel', Xticks);
%     PostFigure([0.8, 2.2, 0, 1], [], 'Mean P(m)', 'Color Wheel');
%     subplot(2,2,2);
%     plot(1:2, meanPm(:,2));
%     set(gca,'XtickLabel', Xticks);
%     PostFigure([0.8,2.2, 0, 1], [], 'Mean P(m)', 'Grey Wheel');
%     subplot(2,2,3);
%     plot(1:2, meanCWA(:,1)); % color-wheel
%     set(gca,'XtickLabel', Xticks);
%     PostFigure([0.8,2.2, -0.03, 0.5], [], 'Mean P(CW)', 'Color Wheel');
%     subplot(2,2,4);
%     plot(1:2, meanCWA(:,2));
%     set(gca,'XtickLabel', Xticks);
%     PostFigure([0.8,2.2, -0.03, 0.5], [], 'Mean P(CW)', 'Grey Wheel');
% 
% end

% 
% if fitIMSim
%     disp('      X       Y          s       kappa      kappaf      r');
%     disp(mean(IMSimparms, 1));
%     disp('      X       Y          s       kappa      kappaf      r');
%     disp(std(IMSimparms, 1));
% end

%%% Save results
if E.saveResults == 1
    fid = fopen('IMSim.RetroCueDecayInter.dat', 'w');
    for id = 1:E.nsubj
        for ri = 1:length(RI)
            for cueing = 1:2
                fprintf(fid, '%d %d %d %d  ', id, RI(ri), cueing, Mdevobs(id, cueing, wheel));
                if (fitMM == 1), fprintf(fid, '%d %d %d %d', MMtranspos(id, cueing, ri), MMguessing(id, cueing, ri), MMcwattraction(id, cueing, ri), MMSD(id, cueing, ri)); end
                fprintf(fid, '\n');
            end
        end
    end
    fclose(fid);

end



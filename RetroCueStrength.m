function [] = RetroCueStrength(Model, fitMM)
% Simulation of retro-cue effects in continuous reproduction, varying strengthening
% Retro-cue effect even at strength = 0. The valid-cue benefit is probaly
% due to a head start of retrieval: It goes away when CSI = 0.01;
% the invalid-cue cost is due to erasure of the high-precision FoA in 1/N trials.
% with setsize 1, there is a cue benefit at small (or 0) strengthening when
% CSI = 1; it goes away when CSI = 0.01.

% Problems
% (1) Strengthening per se has no beneficial effect - it only impairs memory for the not-cued items
% (2) Refreshing pulls FoA from last-encoded element -> cost relative to condition with no refreshing cues


global P
global E
global C

E.PreRetro = 2;
E.cuevalidity = 2/3;
P.cwinter = 0; 
setsize = 6;
Strength = [0, 0.5, 1, 2, 3, 4];
option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

% initialize parameters of mixture model, IM, observed mean Deviation

Mdevobs = NaN(E.nsubj, 3, length(Strength));  % mean observed deviation of responses from target feature for subjects, cueing conditions, and strengthening
Mrt = NaN(E.nsubj, 3, length(Strength));      % mean RT
MMSD = NaN(E.nsubj, 3, length(Strength));     % Mixture Model SD parameter
MMguessing = NaN(E.nsubj, 3, length(Strength));  % Mixture Model Guessing parameter
MMtranspos = NaN(E.nsubj, 3, length(Strength));  % Mixture Model Transposition parameter (swap error proportion)
MMcwattraction = NaN(E.nsubj, 3, length(Strength));  % Mixture Model colorwheel-attraction strength parameter
Mwact = NaN(E.nsubj, 3, length(Strength));       % mean activation in the binding weight matrix

% pre-define container matrices for memory arrays, target features, and
% responses
Array = zeros(E.nsubj*3*length(Strength)*E.ntrials, setsize);
Target = zeros(1,E.nsubj*3*length(Strength)*E.ntrials);
Response = zeros(1,E.nsubj*3*length(Strength)*E.ntrials);

% for each set-size level, for each trial, generate a vector of 360 color
% values coding the colors on the wheel
[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x setsize x [1:360]

tcount = 1; %trial count

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);


    for str = 1:length(Strength)

        P.cueingStrength = Strength(str);

        % Initialize container vectors
        fdistance = zeros(1,E.ntrials);      % feature distance between response and target
        rt = zeros(1,E.ntrials);             % response time
        wact = zeros(1,E.ntrials);           % activation (summed strength of bindings) of binding weight matrix
        Probedpos = zeros(E.ntrials,1);      % Number of the tested (probed) spatial position
        Pangle = zeros(E.ntrials,setsize);   % spatial angles of item positions in the array
        Cangle = zeros(E.ntrials,setsize+1); % color angles in the color wheel
        Targ = zeros(E.ntrials,1);           % Target
        Resp = zeros(E.ntrials,1);           % Response
        Setsize = zeros(E.ntrials,1);        % Set size

        for cueing = 1:3  % neutral, valid, invalid

            for trial = 1:E.ntrials

                output = Model(P, setsize, cueing);   % here the model is run!

                Array(tcount,:) = [output.F(1,1:setsize), zeros(1, setsize-setsize)];  % record of the array on this trial (fill up the rest until max(setsize) with zeros
                Target(tcount) = output.F(1,1);         % target feature
                Response(tcount) = output.response;     % response feature
                rt(trial) = output.rt;                  % response time
                fdistance(trial) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
                wact(trial) = sum(sum(output.wx));      % sum of activation in weight matrix -> CDA?
                tcount = tcount+1;                      % trial counter is incremented

                %collect data for further modeling with Mixture Model or IM
                Probedpos(trial) = output.L(1);         % probed position
                Pangle(trial,:) = output.L(1:setsize);  % spatial position angles
                Cangle(trial,1:setsize) = output.F(1,1:setsize);  % array feature angles
                Cangle(trial,setsize+1) = output.CWcolor;  % color-wheel feature angles
                Targ(trial) = output.F(1,1);            % target feature
                Resp(trial) = output.response;          % response feature
                Setsize(trial) = setsize;               % set size

            end

            Mdevobs(id, cueing, str) = mean(abs(fdistance));  %mean deviation (averaged over trials)
            Mrt(id, cueing, str) = mean(rt);
            Mwact(id, cueing, str) = mean(wact);

            ssData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting
            ssData.cueing = repmat(cueing, length(Response), 1);
            ssData.preretro = E.PreRetro;
            ssD.setsize = Setsize;
            ssD.response = Resp;
            ssD.L = round(C.Location(Pangle));
            ssD.Color = Cangle;
            ssD.cueing = repmat(cueing, length(Response), 1);
            ssD.preretro = E.PreRetro;

            if str == 1
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
                disp('    id        strength  cueing    Deviance  SD(deg)   P(guess)  P(trans)  P(CW)');
                disp([id, Strength(str), cueing, MMloglik/1000, MMparms]);
                MMSD(id, cueing, str) = MMparms(1);
                MMguessing(id, cueing, str) = MMparms(2);
                if npar > 2, MMtranspos(id, cueing, str) = MMparms(3); end
                if npar > 3, MMcwattraction(id, cueing, str) = MMparms(4); end
            end

            disp('    id   strength cueing ');
            disp([id, str, cueing]);

        end % for cueing
    end %for str

end  % for ID

% Plot Mean(Deviation) as functions of set size

legendtext = {'Neutral', 'Valid', 'Invalid'};

% Stand-Alone Figure of Errors
PreFigure;
plotvector = squeeze(mean(Mdevobs,1))';
plot(Strength, plotvector);
PostFigure([0, max(Strength)+0.02, 0, 80], 'Strength', 'Deviation (Deg)', [], legendtext);

% Errors together with Mixture-Model parameters and CDA simulation
PreFigure;
subplot(2,2,1);
Mdevobs(:,3,1) = NaN;  % invalid cue at set size 1
plotvector = squeeze(mean(Mdevobs,1))';
plot(Strength, plotvector);
PostFigure([-0.03 max(Strength)+0.02, 0, 1.05*max(max(plotvector))], 'Strength', 'Deviation', 'Mean Deviation', legendtext);
subplot(2,2,2);
Mrt(:,3,1) = NaN;
plotvector = squeeze(mean(Mrt,1))';
plot(Strength, plotvector);
PostFigure([-0.03 max(Strength)+0.02, 0, 1.05*max(max(plotvector))], 'Strength', 'RT(ms)', 'Mean RT', legendtext);
subplot(2,2,3);
Mwact(:,3,1) = NaN;
plotvector = squeeze(mean(Mwact,1))';
plot(Strength, plotvector);
PostFigure([-0.03 max(Strength)+0.02, 0, 1.1*max(max(plotvector))], 'Strength', 'W.act', 'Weight Activity', legendtext);

% Plot response distributions
meanDeviation = ResponseDistrib(Array, Target, Response);
disp(meanDeviation);


% Plot Mixture Model Parameters over Setsize
if fitMM
    MMPm = 1 - MMtranspos - MMguessing - MMcwattraction;
    meanMMPm = squeeze(mean(MMPm,1));
    meanK = meanMMPm * setsize;
    PreFigure
    subplot(3,2,1);
    plot(Strength, squeeze(mean(MMSD,1))');
    PostFigure([-0.03 max(Strength)+0.02, 0, max(max(mean(MMSD,1)))+0.5], 'Strength', 'Mean SD', 'SD from Bays Mixture', {'Neutral', 'Valid', 'Invalid'});
    subplot(3,2,2);
    plot(Strength, meanK');
    PostFigure([-0.03 max(Strength)+0.02, 0, 6], 'Strength', 'K', 'K from Bays Mixture', {'Neutral', 'Valid', 'Invalid'});
    subplot(3,2,3);
    plot(Strength, squeeze(mean(MMPm,1))');
    PostFigure([-0.03 max(Strength)+0.02, 0, 1], 'Strength', 'Mean P(m)', 'P(mem)');
    subplot(3,2,4);
    plot(Strength, squeeze(mean(MMguessing,1))');
    PostFigure([-0.03 max(Strength)+0.02, 0, 0.5], 'Strength', 'Mean P(guess)', 'P(guess)');
    subplot(3,2,5);
    plot(Strength, squeeze(mean(MMtranspos,1))');
    PostFigure([-0.03 max(Strength)+0.02, 0, 1], 'Strength', 'Mean P(trans)', 'P(transpos)');
    subplot(3,2,6);
    plot(Strength, squeeze(mean(MMcwattraction,1))');
    PostFigure([-0.03 max(Strength)+0.02, 0, 1], 'Strength', 'Mean P(wheel)', 'P(wheel attraction)');
end

%%% Save results
if E.saveResults == 1
    fid = fopen('IMSim.RetroCueStrength.dat', 'w');
    for id = 1:E.nsubj
        for str = 1:length(Strength)
            for cueing = 1:3
                fprintf(fid, '%d %d %d %d  ', id, Strength(str), cueing, Mdevobs(id, cueing, str));
                for ii = 1:length(indVar)
                    fprintf(fid, '%d ', ParX(id, ii));
                end
                if (fitMM == 1), fprintf(fid, '%d %d %d %d', MMtranspos(id, cueing, str), MMguessing(id, cueing, str), MMcwattraction(id, cueing, str), MMSD(id, cueing, str)); end
                fprintf(fid, '\n');
            end
        end
    end
    fclose(fid);
end




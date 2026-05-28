function [] = CueTargetInterval(Model, setsize, fitMM, fitIMSim)
% Simulation of retro-cue in continuous reproduction, varying cue-target
% interval (Souza et al., 2016)

global P
global E
global C

E.PreRetro = 2;
E.cuevalidity = 1;
option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);
fitIM = 0;
CTI = [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]; 

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% initialize parameters of mixture model, IM, observed mean Deviation

Mdevobs = NaN(E.nsubj, 3, length(CTI));  % mean observed deviation of responses from target feature for subjects, cueing conditions, and set size
Mrt = NaN(E.nsubj, 3, length(CTI));      % mean RT
MMSD = NaN(E.nsubj, 3, length(CTI));     % Mixture Model SD parameter
MMguessing = NaN(E.nsubj, 3, length(CTI));  % Mixture Model Guessing parameter
MMtranspos = NaN(E.nsubj, 3, length(CTI));  % Mixture Model Transposition parameter (swap error proportion)
MMcwattraction = NaN(E.nsubj, 3, length(CTI));  % Mixture Model colorwheel-attraction strength parameter
Mwact = NaN(E.nsubj, 3, length(CTI));       % mean activation in the binding weight matrix

IMparms = zeros(E.nsubj, 6);     % Parameters for Interference Mdel
IMSimparms = zeros(E.nsubj, 6);   % Simulation parameters for IM

% pre-define container matrices for memory arrays, target features, and
% responses
Array = zeros(E.nsubj*3*length(CTI)*E.ntrials, setsize);
Target = zeros(1,E.nsubj*3*length(CTI)*E.ntrials);
Response = zeros(1,E.nsubj*3*length(CTI)*E.ntrials);

% for each array item, for each trial, generate a vector of 360 color
% values coding the colors on the wheel
[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x setsize x [1:360]


% generate parameters with individual differences
ParX = CreateIndDiff;

tcount = 1; %trial count

for id = 1:E.nsubj
    
        
    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);
    
    for cti = 1:length(CTI)
        
        E.CTI(2:3) = CTI(cti);  % overwrite default values
        E.CQI(2:3) = CTI(cti);
        
        % Initialize container vectors
        fdistance = zeros(1,E.ntrials);  % feature distance between response and target
        rt = zeros(1,E.ntrials);         % response time
        wact = zeros(1,E.ntrials);       % activation (summed strength of bindings) of binding weight matrix
        Probedpos = zeros(E.ntrials,1);  % Number of the tested (probed) spatial position
        Pangle = zeros(E.ntrials,setsize);  % spatial angles of item positions in the array
        Cangle = zeros(E.ntrials,setsize+1); % color angles in the color wheel
        Targ = zeros(E.ntrials,1);                % Target
        Resp = zeros(E.ntrials,1);                % Response
        Setsize = zeros(E.ntrials,1);             % Set size
        
        for cueing = 1:2  % neutral, valid - there were no invalid cues in the experiment of Souza et al. (2016)
            
            for trial = 1:E.ntrials
                
                output = Model(P, setsize, cueing);   % here the model is run!
                
                Array(tcount,:) = output.F(1:setsize);  % record of the array on this trial (fill up the rest until max(setsize) with zeros
                Target(tcount) = output.F(1);    % target feature
                Response(tcount) = output.response;   % response feature
                rt(trial) = output.rt;            % response time
                fdistance(trial) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
                wact(trial) = sum(sum(output.wx)); % sum of activation in weight matrix -> CDA?
                tcount = tcount+1;   % trial counter is incremented
                
                %collect data for further modeling with Mixture Model or IM
                Probedpos(trial) = output.L(1);    % probed position
                Pangle(trial,:) = output.L(1:setsize);  % spatial position angles
                Cangle(trial,1:setsize) = output.F(1:setsize);  % array feature angles
                Cangle(trial,setsize+1) = output.CWcolor;  % color-wheel feature angles
                Targ(trial) = output.F(1);           % target feature
                Resp(trial) = output.response;       % response feature
                Setsize(trial) = setsize;            % set size
                
            end
            
            Mdevobs(id, cueing, cti) = mean(abs(fdistance));  %mean deviation (averaged over trials)
            Mrt(id, cueing, cti) = mean(rt);
            Mwact(id, cueing, cti) = mean(wact);
            
            ssData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting
            ssData.cueing = repmat(cueing, length(Response), 1);
            ssData.preretro = E.PreRetro;
            ssD.setsize = Setsize;
            ssD.response = Resp;
            ssD.L = round(C.Location(Pangle));
            ssD.Color = Cangle;
            ssD.cueing = repmat(cueing, length(Response), 1);
            ssD.preretro = E.PreRetro;
            
            if cti == 1
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
                disp([id, CTI(cti), itercount, MMloglik]);
                MMSD(id, cueing, cti) = MMparms(1);
                MMguessing(id, cueing, cti) = MMparms(2);
                if npar > 2, MMtranspos(id, cueing, cti) = MMparms(3); end
                if npar > 3, MMcwattraction(id, cueing, cti) = MMparms(4); end
            end
            
        end % for cueing
    end %for setsize
    
    %fit IM (to all set size conditions jointly)
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
    
    %fit IMFit
    if fitIMSim
        startparms = [0.05, 1.5, 3, 10, 20, 0.5];  %X, Y, s, P.kappa, P.kappafocus, Creduction
        npar = 6;
        lb = zeros(1,npar); ub = [5, 5, 20, 90, 90, 1];
        IMSimloglik = 500000;
        itercount = 0;
        while IMSimloglik > 400000
            [IMSimparms(id,:), IMSimloglik] = fminsearchbnd(@(x) IMFit(x, D, 2), startparms, lb, ub, option);
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

% Plot Mean(Deviation) as functions of set size

legendtext = {'Neutral', 'Valid', 'Invalid'};

% Stand-Alone Figure of Errors
PreFigure;
plotvector = squeeze(mean(Mdevobs,1))';
plot(CTI, plotvector);
PostFigure([0, max(CTI)+0.1, 0, 1.05*max(max(plotvector))], 'Cue Target Interval', 'Deviation (Deg)', [], legendtext);

% Plot response distributions
%meanDeviation = ResponseDistrib(Array, Target, Response);
%disp(meanDeviation);

% Plot Mixture Model Parameters over Setsize
if fitMM
    MMPm = 1 - MMtranspos - MMguessing - MMcwattraction;
    meanMMPm = squeeze(mean(MMPm,1));
    meanK = bsxfun(@times, meanMMPm, 1:8);
    PreFigure;
    subplot(3,2,1);
    plot(CTI, squeeze(mean(MMSD,1))');
    PostFigure([0, max(CTI)+0.1, 0, max(max(mean(MMSD,1)))+0.5], 'Cue Target Interval', 'Mean SD', 'SD from Bays Mixture', {'Neutral', 'Valid', 'Invalid'});
    subplot(3,2,2);
    plot(CTI, meanK');
    PostFigure([0, max(CTI)+0.1, 0, 6], 'Cue Target Interval', 'K', 'K from Bays Mixture', {'Neutral', 'Valid', 'Invalid'});
    subplot(3,2,3);
    plotCTI, (squeeze(mean(MMPm,1))');
    PostFigure([0, max(CTI)+0.1, 0, 1], 'Cue Target Interval', 'Mean P(m)', 'P(mem)');
    subplot(3,2,4);
    plot(CTI, squeeze(mean(MMguessing,1))');
    PostFigure([0, max(CTI)+0.1, 0, 0.5], 'Cue Target Interval', 'Mean P(guess)', 'P(guess)');
    subplot(3,2,5);
    plot(CTI, squeeze(mean(MMtranspos,1))');
    PostFigure([0, max(CTI)+0.1, 0, 1], 'Cue Target Interval', 'Mean P(trans)', 'P(transpos)');
    subplot(3,2,6);
    plot(CTI, squeeze(mean(MMcwattraction,1))');
    PostFigure([0, max(CTI)+0.1, 0, 1], 'Cue Target Interval', 'Mean P(wheel)', 'P(wheel attraction)');
end

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
    fid = fopen('IMSim.CTI.dat', 'w');
    for id = 1:E.nsubj
        for cti = 1:length(CTI)
            for cueing = 1:3
                fprintf(fid, '%d %d %d %d  ', id, cti, cueing, Mdevobs(id, cueing, cti));
                if (fitMM == 1), fprintf(fid, '%d %d %d %d', MMtranspos(id, cueing, cti), MMguessing(id, cueing, cti), MMcwattraction(id, cueing, cti), MMSD(id, cueing, cti)); end
                fprintf(fid, '\n');
            end
        end
    end
    fclose(fid);
end




function [] = ParameterSensitivity2(Model, ParName, Values)
% Simulation of effect of varying a parameter in a simultaneous-presentation
% continuous-reproduction paradigm with a single test

global P
global E
global C

setsize = 6;
E.presentation = 1;  % sequential!
E.outsize = 1;
fitMM = 1;
fitSDM = 0;
option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end


[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x E.maxsetsize x [1:360]

tcount = 1; %trial count

% initialize parameters of mixture model, IM, observed mean Deviation

Mdevobs = NaN(E.nsubj, length(Values));  % mean observed deviation of responses from target feature for subjects, cueing conditions, and set size
Mrt = NaN(E.nsubj,  length(Values));      % mean RT
MMSD = NaN(E.nsubj,  length(Values));     % Mixture Model SD parameter
MMguessing = NaN(E.nsubj,  length(Values));  % Mixture Model Guessing parameter
MMtranspos = NaN(E.nsubj,  length(Values));  % Mixture Model Transposition parameter (swap error proportion)
MMcwattraction = NaN(E.nsubj,  length(Values));  % Mixture Model colorwheel-attraction strength parameter
SDMc = NaN(E.nsubj,  length(Values));        % SDM parameter c
SDMkappa = NaN(E.nsubj,  length(Values));    % SDM parameter Kappa
SDMa = NaN(E.nsubj,  length(Values));        % SDM parameter a
SDMs = NaN(E.nsubj,  length(Values));        % SDM parameter s
Mwact = NaN(E.nsubj,  length(Values));       % mean activation in the binding weight matrix

% pre-define container matrices for memory arrays, target features, and
% responses
Array = zeros(E.nsubj*length(Values)*E.ntrials, setsize);
Target = zeros(1,E.nsubj*length(Values)*E.ntrials);
Response = zeros(1,E.nsubj*length(Values)*E.ntrials);

% generate parameters with individual differences

PP = P; % keep original Parameter structure

for parIdx = 1:length(Values)

    P = PP; 
    eval(['P.', ParName, ' = Values(parIdx);']);
    ParX = CreateIndDiff;

    for id = 1:E.nsubj

        % extract parameter values for each subject - for those parameters that vary between subjects
        for ii = 1:length(C.indVar)
            eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
        end

        % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
        CreateStimuli;
        CreateMapping(E.calibrateAmp==2);

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

        for trial = 1:E.ntrials

            output = Model(P, setsize, 1);   % here the model is run!

            Array(tcount,:) = output.F(1,1:setsize);  % record of the array on this trial (fill up the rest until max(setsize) with zeros
            Target(tcount) = output.F(1,1);    % target feature
            Response(tcount) = output.response;   % response feature
            rt(trial) = output.rt;            % response time
            fdistance(trial) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
            wact(trial) = sum(sum(output.wx)); % sum of activation in weight matrix -> CDA?
            tcount = tcount+1;   % trial counter is incremented

            %collect data for further modeling with Mixture Model or IM
            Probedpos(trial) = output.L(1);    % probed position
            Pangle(trial,:) = output.L(1:setsize);  % spatial position angles
            Cangle(trial,1:setsize) = output.F(1,1:setsize);  % array feature angles
            Cangle(trial,setsize+1) = output.CWcolor;  % color-wheel feature angles
            Targ(trial) = output.F(1,1);           % target feature
            Resp(trial) = output.response;       % response feature
            Setsize(trial) = setsize;            % set size

        end

        Mdevobs(id, parIdx) = nanmean(abs(fdistance));  %mean deviation
        Mrt(id, parIdx) = mean(rt);
        Mwact(id, parIdx) = mean(wact);

        Data = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting

        % fit Mixture Model
        if fitMM
            startparms = [15, .1, .1, .1];
            lb = [eps, 0, 0, 0]; ub = [90, 1, 1, 1];
            npar = 4;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture, 4 -> Souza & Oberauer mixture (including color-wheel attraction)
            MMloglik = 500000;
            itercount = 0;
            while MMloglik > 400000
                [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, Data, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                itercount = itercount + 1;
            end
            disp('    ID        parameter iteration Deviance SD(deg)   P(guess)  P(trans)   P(CW)');
            disp([id, Values(parIdx), itercount, MMloglik/1000, MMparms]);
            MMSD(id, parIdx) = MMparms(1);
            MMguessing(id, parIdx) = MMparms(2);
            if npar > 2, MMtranspos(id, parIdx) = MMparms(3); end
            if npar > 3, MMcwattraction(id, parIdx) = MMparms(4); end
        end

        if fitSDM
            Data.feature = pi*Cangle./180;
            startparms = [4, 5, 1, 2];  % C, kappa, A, s
            lb = [eps, eps, 0, 0]; ub = [30, 70, 10, 10];
            SDMloglik = 500000;
            itercount = 0;
            while SDMloglik > 400000
                [SDMparms, SDMloglik] = fminsearchbnd(@(x) SDM(x, Data, 2), startparms, lb, ub, option);
                itercount = itercount + 1;
            end
            disp('    ID        parameter iteration Deviance  C        Kappa      A         s');
            disp([id, Values(parIdx), itercount, SDMloglik/1000, SDMparms]);
            SDMc(id, parIdx) = SDMparms(1);
            SDMkappa(id, parIdx) = SDMparms(2);
            SDMa(id, parIdx) = SDMparms(3);
            SDMs(id, parIdx) = SDMparms(4);

        end

        if fitMM == 0 && fitSDM == 0
            disp('ID  Parameter  Value');
            disp([mat2str(id), '   ' ParName, '   ', mat2str(Values(parIdx))]);
        end

    end %for ID
end  % for parIdx

% Plot Mean(Deviation) as function of set size

Dash = find(ParName=='_');
if ~isempty(Dash)
    NewParName = [ParName(1:(Dash-1))];
    for letter = (Dash+1):length(ParName)
        NewParName = [NewParName, '_', ParName(letter)]; 
    end
    ParName = NewParName;
end


% Stand-alone figures of errors
legendtext = {'In=1', 'In=2', 'In=3','In=4', 'In=5', 'In=6', 'In=7', 'In=8'};
PreFigure;
plotvector = mean(Mdevobs);  
plot(Values, plotvector);
PostFigure([min(Values)-0.05*max(Values), 1.05*max(Values), 0, 1.05*max(max(plotvector))], ParName, 'Deviation (Deg)', [], legendtext);


% Plot Mixture Model Parameters over Parameter Values
if fitMM
    MMPm = 1 - MMtranspos - MMguessing - MMcwattraction;
    meanMMPm = squeeze(mean(MMPm,1));
    K = meanMMPm .* setsize;
    PreFigure;
    subplot(3,2,1);
    plot(Values, squeeze(mean(MMSD,1))');
    PostFigure([min(Values)-0.2*abs(min(Values)), 1.2*max(Values), 0, max(max(mean(MMSD,1)))+0.5], ParName, 'Mean SD', 'SD from Bays Mixture');
    subplot(3,2,2);
    plot(Values, K');
    PostFigure([min(Values)-0.2*abs(min(Values)), 1.2*max(Values), 0, 6], ParName, 'K', 'K from Bays Mixture');
    subplot(3,2,3);
    plot(Values, squeeze(mean(MMPm,1))');
    PostFigure([min(Values)-0.2*abs(min(Values)), 1.2*max(Values), 0, 1], ParName, 'Mean P(m)', 'P(mem)');
    subplot(3,2,4);
    plot(Values, squeeze(mean(MMguessing,1))');
    PostFigure([min(Values)-0.2*abs(min(Values)), 1.2*max(Values), 0, 1], ParName, 'Mean P(guess)', 'P(guess)');
    subplot(3,2,5);
    plot(Values, squeeze(mean(MMtranspos,1))');
    PostFigure([min(Values)-0.2*abs(min(Values)), 1.2*max(Values), 0, 0.5], ParName, 'Mean P(trans)', 'P(transpos)');
    subplot(3,2,6);
    plot(Values, squeeze(mean(MMcwattraction,1))');
    PostFigure([min(Values)-0.2*abs(min(Values)), 1.2*max(Values), 0, 1], ParName, 'Mean P(wheel)', 'P(wheel attraction)');
end


% Plot SDM Parameters over Parameter Values
if fitSDM
    PreFigure;
    subplot(2,2,1);
    plot(Values, squeeze(mean(SDMc,1))');
    PostFigure([min(Values)-0.2*abs(min(Values)), 1.2*max(Values), 0, max(max(mean(SDMc,1)))+1], ParName, 'Mean C', 'C from SDM');
    subplot(2,2,2);
    plot(Values, squeeze(mean(SDMkappa,1))');
    PostFigure([min(Values)-0.2*abs(min(Values)), 1.2*max(Values), 0, max(max(mean(SDMkappa,1)))+1], ParName, 'Mean Kappa', 'Kappa from SDM');
    subplot(2,2,3);
    plot(Values, squeeze(mean(SDMa,1))');
    PostFigure([min(Values)-0.2*abs(min(Values)), 1.2*max(Values), 0, max(max(mean(SDMa,1)))+0.5], ParName, 'Mean A', 'A from SDM');
    subplot(2,2,4);
    plot(Values, squeeze(mean(SDMs,1))');
    PostFigure([min(Values)-0.2*abs(min(Values)), 1.2*max(Values), 0, max(max(mean(SDMs,1)))+1], ParName, 'Mean s', 's from SDM');
end



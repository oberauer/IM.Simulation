function [] = ParameterSensitivityCD(Model, ParName, Values)
% Simulation of effect of varying a parameter in a simultaneous-presentation
% change-detection paradigm with a single test

global P
global E
global C

E.presentation = 1;  
E.outsize = 1;
E.test = 2; 
E.wheel = 0; 
E.material = 2; 
setsize = 3;
cueing = 1; 
nTrials = round(E.ntrials./4);

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

IMprepareRecog; % set up criterion for expected size of change 

% Initializing some container matrices
Pyes = zeros(E.nsubj, 3, length(Values));  % Probability of saying "Yes" (="Same") for each subject, probe type, level of parameter
PC = zeros(E.nsubj, 3, length(Values));    % Proportion correct
RT = zeros(E.nsubj, 3, length(Values));    % Response time

Ptype = [1, 1, 2, 3];   % levels of the probetype variable 2 x positive, 1 x new, 1 x intrusion
Design = fullfact(4);   % probetype
nCells = size(Design, 1);      % number of design cells (crossing cueing, probetype array 1, and probetype array 2)

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

        % Prepare matrices for results
        Conditionvector = repmat(1:nCells, nTrials);
        Conditionvector = Conditionvector(randperm(length(Conditionvector)));  %shuffle condition vector
        ConditionCount = zeros(1,nCells);
        response = zeros(nTrials, nCells);
        rt = zeros(nTrials, nCells);

        % running counter of trials in each condition
        for trial = 1:(nTrials*nCells)
            condition = Conditionvector(trial);         % pick the condition of this trial
            ConditionCount(condition) = ConditionCount(condition) + 1;  % increment trial count for the current trial's condition
            E.ptype = Ptype(Design(condition, 1));             % determine the probetype from the design matrix
            output = Model(P, setsize, cueing);   % run model on 1 trial, returns predictions (output is a structure with lots of variables in it)
            response(ConditionCount(condition), condition) = output.response(1,:);  % the first entry of response is the actual response
            %disp([cueing, E.ptype, output.response(1,:)]);
            rt(ConditionCount(condition), condition) = output.rt;    % response time
        end

        % now loop over the 4 design cells to read out the summary
        % statistics of simulated data in each cell
        for condition = 1:4
            ptype = Ptype(Design(condition,1));
            Pyes(id, ptype, parIdx) = mean(2-response(:,condition));  % Yes/No: response = 1/2
            if (ptype == 1)
                PC(id, ptype, parIdx) = Pyes(id, ptype, parIdx);
            else
                PC(id, ptype, parIdx) = 1-Pyes(id, ptype, parIdx);
            end
            RT(id, ptype, parIdx) = mean(rt(:,condition))./1000;
        end

        disp('ID  Parameter  Value    PC');
        disp([mat2str(id), '   ' ParName, '   ', mat2str(Values(parIdx)), '   ', mat2str(round(mean(PC(id,:,parIdx)),2))]);

    end %for ID
end  % for parIdx

% Plot Mean(PC) as function of set size

Dash = find(ParName=='_');
if ~isempty(Dash)
    NewParName = [ParName(1:(Dash-1))];
    for letter = (Dash+1):length(ParName)
        NewParName = [NewParName, '_', ParName(letter)]; 
    end
    ParName = NewParName;
end


disp(squeeze(mean(PC)));

K = setsize.*((squeeze(mean(Pyes(:,1,:), 1))-squeeze(mean(Pyes(:,2,:), 1)))');

legendtext = {'Positive', 'New', 'Intrusion'}; 
PreFigure;
subplot(1,3,1);
plotvector = squeeze(mean(PC));
plot(Values, plotvector);
PostFigure([min(Values)-0.05*max(Values), 1.05*max(Values), 0, 1], ParName, 'P(correct)', [], legendtext);
subplot(1,3,2);
plotvector = squeeze(mean(RT));
plot(Values, plotvector);
PostFigure([min(Values)-0.05*max(Values), 1.05*max(Values), 0, 1.05*max(max(plotvector))], ParName, 'RT(s)', [], legendtext);
subplot(1,3,3);
plot(Values, K);
PostFigure([min(Values)-0.05*max(Values), 1.05*max(Values), 0, 1.05*max(max(K))], ParName, 'K');

output = 0; 

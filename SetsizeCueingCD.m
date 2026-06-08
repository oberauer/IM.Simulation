function D = SetsizeCueingCD(Model)
% Simulation of Interaction of retro-cueing with set-size in Change Detection
% Experiment 1 of Souza & Oberauer (2014,JEP:HPP), varying set size (1 to 6) and
% retro-cueing (no cue vs. valid cue). 

global P
global E
global C

E.test = 2;      % change detection
E.PreRetro = 2;  % this is all retro-cue
E.cuevalidity = 1;
E.maxsetsize = 6; 
nTrials = round(E.ntrials./4);

C.nstim = 9;    % 9 highly discriminable stimuli were used in the experiment
IMprepareRecog; % set up criterion for expected size of change 

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

% Initializing some container matrices
Pyes = zeros(E.nsubj, 3, E.maxsetsize, 2);  % Probability of saying "Yes" (="Same") for each subject, probe type, setsize, and cueing condition
PC = zeros(E.nsubj, 3, E.maxsetsize, 2);    % Proportion correct
RT = zeros(E.nsubj, 3, E.maxsetsize, 2);    % Response time

Ptype = [1, 1, 2, 3];   % levels of the probetype variable 2 x positive, 1 x new, 1 x intrusion
Design = fullfact([4, 2]);   % probetype x cueing
nCells = size(Design, 1);  % number of design cells (crossing cueing with probetype)

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);

    for setsize = 1:E.maxsetsize

        C.nloc = setsize;  % items were dispersed maximally in space

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
            if setsize==1, E.ptype = min(2, E.ptype); end      % no intrusion probes in set size 1
            cueing = Design(condition, 2);
            output = Model(P, setsize, cueing);   % run model on 1 trial, returns predictions (output is a structure with lots of variables in it)
            response(ConditionCount(condition), condition) = output.response(1,:);  % the first entry of response is the actual response
            %disp([cueing, E.ptype, output.response(1,:)]);
            rt(ConditionCount(condition), condition) = output.rt;    % response time
        end

        % now loop over the 8 design cells to read out the summary
        % statistics of simulated data in each cell
        for condition = 1:size(Design,1)
            ptype = Ptype(Design(condition,1));
            cueing = Design(condition,2);
            Pyes(id, ptype, setsize, cueing) = mean(2-response(:,condition));  % Yes/No: response = 1/2
            if (ptype == 1)
                PC(id, ptype, setsize, cueing) = PC(id, ptype, setsize, cueing) + Pyes(id, ptype, setsize, cueing)./2;  % divide by 2 because there are 2 conditions for positive probes
            else
                PC(id, ptype, setsize, cueing) = 1-Pyes(id, ptype, setsize, cueing);
            end
            RT(id, ptype, setsize, cueing) = mean(rt(:,condition));
        end

        disp(['      ID      Setsize   PC no-cue PC cue    RT']);
        disp([id, setsize, mean(PC(id, :, setsize, 1)), mean(PC(id, :, setsize, 2)), mean(mean(RT(id, :, setsize, :)))]);

    end

end  % for ID

disp(squeeze(mean(PC(:,:,:,1)))); % K by setsize and cueing condition

K = (1:E.maxsetsize).*((squeeze(mean(Pyes(:,1,:,:), 1))-squeeze(mean(Pyes(:,2,:,:), 1)))'); % P(yes|pos) - P(yes|new)

% Accuracy and RT by set size and cueing condition
PreFigure;
subplot(2,2,1);
plotvector = (2*squeeze(mean(PC(:,1,:,:),1)) + squeeze(mean(PC(:,2,:,:),1)) + squeeze(mean(PC(:,3,:,:),1)))./4;  % average over id and ptype, 2x weight for positive probes
plot(1:E.maxsetsize, plotvector);  % average over id and ptype
PostFigure([0.5, E.maxsetsize + 0.5, 0, 1], 'Set Size', 'P(correct)', 'Accuracy, All Probe Types', {'No Cue', 'Cue'});
subplot(2,2,2);
plot(1:E.maxsetsize, squeeze(mean(mean(PC(:,1:2,:,:),2),1)));  % average over id and ptype, excluding intrusions
PostFigure([0.5, E.maxsetsize + 0.5, 0, 1], 'Set Size', 'P(correct)', 'Accuracy, Pos. and New', {'No Cue', 'Cue'});
subplot(2,2,3);
plotvector = (2*squeeze(mean(RT(:,1,:,:))) + squeeze(mean(RT(:,2,:,:))) + squeeze(mean(RT(:,3,:,:))))./4;
plot(1:E.maxsetsize, plotvector);  % average over id and ptype
PostFigure([0.5, E.maxsetsize + 0.5, 0, 1.2*max(plotvector(:))], 'Set Size', 'RT(s)', 'RT, All Probe Types', {'No Cue', 'Cue'});
subplot(2,2,4);
plot(1:E.maxsetsize, squeeze(mean(K)));
PostFigure([0, E.maxsetsize + 0.5, 0.5, 4], 'Set Size', 'K', 'Cowan K', {'No Cue', 'Cue'});

% Accuracy and RT by set size and probe type, no-cue condition only
PreFigure;
subplot(1,2,1);
plotvector = squeeze(mean(PC(:,:,:,1)));  % no-cue condition only
plotvector(3,1) = NaN; % intrusion probe at set size 1
plot(1:E.maxsetsize, plotvector);  
PostFigure([0.5, E.maxsetsize + 0.5, 0, 1], 'Set Size', 'P(correct)', 'Accuracy', {'Positive', 'New', 'Intrusion'});
subplot(1,2,2);
plotvector = squeeze(mean(RT(:,:,:,1)));
plotvector(3,1) = NaN; % intrusion probe at set size 1
plot(1:E.maxsetsize, plotvector);  % no-cue condition only
PostFigure([0.5, E.maxsetsize + 0.5, 0, 1.2*max(plotvector(:))], 'Set Size', 'RT(s)', 'Response Time', {'Positive', 'New', 'Intrusion'});

% Accuracy, separately for Hits and FAs
PreFigure;
subplot(1,3,1);
plotvector = zeros(E.maxsetsize,2); % set size, cue condition
plotvector(:,1) = squeeze(mean(PC(:,1,:,1), 1));  % positive probes, no-cue condition
plotvector(:,2) = squeeze(mean(PC(:,1,:,2), 1));  % positive probes, cue condition
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0.5, 1], 'Set Size', 'P(hit)', [], {'No Cue', 'Cue'});
subplot(1,3,2);
plotvector(:,1) = squeeze(mean(PC(:,2,:,1), 1));  % new probes, no-cue condition
plotvector(:,2) = squeeze(mean(PC(:,2,:,2), 1));  % new probes, cue condition
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0.5, 1], 'Set Size', 'P(CR New)', [], {'No Cue', 'Cue'});
subplot(1,3,3);
plotvector = NaN(E.maxsetsize, 2);
plotvector(2:E.maxsetsize,1) = squeeze(mean(PC(:,3,2:E.maxsetsize,1), 1));  % intrusion probes, no-cue condition
plotvector(2:E.maxsetsize,2) = squeeze(mean(PC(:,3,2:E.maxsetsize,2), 1));  % intrusion probes, cue condition
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0.5, 1], 'Set Size', 'P(CR Intrus)', [], {'No Cue', 'Cue'});

D.PC = PC;
D.Pyes = Pyes;
D.RT = RT;

output = 0; 




function [] = SensoryMemoryCD(Model)
% Simulation of Pratte & Greene (2023): Change detection with 8-item arrays
% followed by a 100% valid retro cue at various SOAs

global P
global C
global E

E.prestime = 0.05; % arrays were shown just 50 ms
E.PreRetro = 2;  % this is all retro-cue
E.CTI(2) = 0.5;  % cue-probe interval in Pratte & Greene (2023)
E.cuevalidity = 1;
E.maxsetsize = 8;
C.nstim = 10;  % 10 highly discriminable stimuli were used in Pratte & Greene
IMprepareRecog;

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

SOA = [33, 83, 167, 300, 500, 1000]./1000;   % SOAs used in Pratte & Greene
setsize = E.maxsetsize;
cueing = E.PreRetro;
nTrials = round(E.ntrials./2);

% Initializing some container matrices
Pyes = zeros(E.nsubj, 2, length(SOA));  % Probability of saying "Yes" (="Same") for each subject, probe type, and array-cue interval
PC = zeros(E.nsubj, 2, length(SOA));    % Proportion correct
RT = zeros(E.nsubj, 2, length(SOA));    % Response time
Strength = zeros(E.nsubj, 2, length(SOA));    % Strength of Content read out of FX

Ptype = [1, 2];   % levels of the probetype variable 1x positive, 1x new (intrusion trials were not presented in Pratte & Greene)
Design = fullfact(2);   % probetype
nCells = size(Design, 1);      % number of design cells (crossing cueing, probetype array 1, and probetype array 2)

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);

    for soa = 1:length(SOA)

        E.RI = SOA(soa);
        Conditionvector = repmat(1:nCells, nTrials);
        Conditionvector = Conditionvector(randperm(length(Conditionvector)));  %shuffle condition vector
        ConditionCount = zeros(1,nCells);
        response = zeros(nTrials, nCells);
        rt = zeros(nTrials, nCells);
        strength = zeros(nTrials, nCells);

        % running counter of trials in each condition
        for trial = 1:(nTrials*nCells)
            condition = Conditionvector(trial);         % pick the condition of this trial
            ConditionCount(condition) = ConditionCount(condition) + 1;  % increment trial count for the current trial's condition
            E.ptype = Ptype(Design(condition, 1));             % determine the probetype from the design matrix
            output = Model(P, setsize, cueing);   % run model on 1 trial, returns predictions (output is a structure with lots of variables in it)
            response(ConditionCount(condition), condition) = output.response(1,:);  % the first entry of response is the actual response
            %disp([cueing, E.ptype, output.response(1,:)]);
            rt(ConditionCount(condition), condition) = output.rt;    % response time
            strength(ConditionCount(condition), condition) = max(output.Strength);  % the maximal strength is the cued item = the target
        end

        % now loop over the 2 design cells to read out the summary
        % statistics of simulated data in each cell
        for condition = 1:2
            ptype = Ptype(Design(condition,1));
            Pyes(id, ptype, soa) = mean(2-response(:,condition));  % Yes/No: response = 1/2
            if (ptype == 1)
                PC(id, ptype, soa) = Pyes(id, ptype, soa);
            else
                PC(id, ptype, soa) = 1-Pyes(id, ptype, soa);
            end
            RT(id, ptype, soa) = mean(rt(:,condition))./1000;
            Strength(id, ptype, soa) = mean(strength(:, condition)); 
        end

        disp(['ID = ', mat2str(id), '; SOA = ', mat2str(SOA(soa))]);

    end


end  % for ID

%%% Plots

% Proportion correct as a function of set size and cueing condition

PreFigure;
plotvector = zeros(length(SOA), 2); 
plotvector(:,1) = squeeze(mean(PC(:,1,:), 1));  % positive probes
plotvector(:,2) = squeeze(mean(PC(:,2,:), 1));  % new probes
K = setsize.*(squeeze(mean(Pyes(:,1,:), 1))-squeeze(mean(Pyes(:,2,:), 1)));
subplot(1,2,1);
plot(SOA, plotvector);
PostFigure([0, 1, 0, 1], 'SOA', 'P(correct)', [], {'Positive', 'New'});
subplot(1,2,2);
plot(SOA, K);
PostFigure([0, 1, 0, setsize], 'SOA', 'K');

PreFigure;
plotvector = squeeze(mean(mean(Strength)));
plot(SOA, plotvector);
PostFigure([0, 1, 0, max(0.01, 1.1*max(plotvector))], 'SOA', 'Strength of Target');

%%% Save results
if E.saveResults == 1
    filename = ['IMSim.SensoryMemoryCD', mat2str(E.material), '.dat'];
    fid = fopen(filename, 'w');
    for id = 1:E.nsubj
        for soa = 1:length(SOA)
            for ptype = 1:2
                fprintf(fid, '%d %d %d %d %d \n', id, SOA(soa), ptype, PC(id, ptype, soa), RT(id, ptype, soa));
            end
        end
    end
    fclose(fid);
end





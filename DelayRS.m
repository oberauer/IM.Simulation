function [] = DelayRS(Model)
% Simulation of Retro-Cue and Response-Selection Delay for Change Detection
% (Experiment 2 of Souza, Rerko, & Oberauer, 2016, JEP:HPP):
% Conditions: (1) Neutral (location cue, probe color, and response question
% are presented simultaneously), (2) Delay (location cue and probe color
% are presented simultaneusly, and response question presented later), (3)
% Retro-Cue (location cue is presented alone, probe color and response
% question presented later)

% Delay has a beneficial effect only on hits (detecting "same"); retro-cue
% has a beneficial effect only on correct rejections (detecting "change").

% Yes-No accumulation before finishing probe retrieval biases towards "no"
% (change) because it is mostly driven by uniform noise, which is predominantly in
% the "change" area of the feature scale. Delay delays that accumulation,
% thereby reducing the tendency to say "no". Probe interference creates a
% bias towards "yes" (same), which the retro-cue reduces because it enables
% retrieval before probe interference -> advantage for "change" probes.

global P
global E
global C
global M

E.PreRetro = 2;  % this is all retro-cue
E.mask = 3;
E.cuevalidity = 1;
setsize = 6;
IMprepareRecog;

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff; 

Pyes = zeros(E.nsubj, 3, 3);  %Probability of saying "yes" ("same") by subject, probe type (positive, new, intrusion) and cue/delay condition (nneutral, delay, retro-cue)
PC = zeros(E.nsubj, 3, 3);
RT = zeros(E.nsubj, 3, 3);
PyesBig = zeros(E.nsubj, 3, 3);  %Probability of saying "yes" ("same") by subject, probe type (positive, new, intrusion) and cue/delay condition (nneutral, delay, retro-cue)
PCbig = zeros(E.nsubj, 3, 3);
RTbig = zeros(E.nsubj, 3, 3);
binBounds = [1, 15:15:180];                 % Bin boundaries for binning the degrees of change among change trials
binBounds(end) = binBounds(end)+0.5;

for id = 1:E.nsubj
    
        % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']); 
    end
    
    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuliCreateStimuli(E.calibrateAmp==2);
    CreateMapping;
    
    M.kappa(setsize) = 10;        % prior for subjective estimate of memory precision (-> Recognition criterion)
    
    Cueing = [1, 5, 2]; % neutral, delay of response selection, valid retro-cue
    CCond = [ones(1,4), 2*ones(1,4), 3*ones(1,4)];   % cueing condition (consecutive numbers that serve as indices into Cueing)
    Ptype = repmat([1 1 2 3], 1, 3);  % 2 x positive, 1 x new, 1 x intrusion
    Design = [CCond', Ptype'];        % full design with 12 cells; first column codes the cue condition; second column codes the probe type
    response = zeros(E.ntrials, 12);  % initiate response matrix (one column for each design cell, one row for each trial)
    rt = zeros(E.ntrials, 12);        % initiate RT matrix
    Conditionvector = repmat(1:12, 1, E.ntrials);  % vector of condition codes for all trials
    Conditionvector = Conditionvector(randperm(length(Conditionvector)));  % shuffle the order of conditions
    ConditionCount = zeros(1,12);     % counter for trials in each of the 12 design cells
    Bin = zeros(E.ntrials, 12);                  % bins for changes of probes
    
    for trial = 1:(12*E.ntrials)
        condition = Conditionvector(trial);    % find current trial's condition in condition vector
        ConditionCount(condition) = ConditionCount(condition) + 1;  % increment trial count for the current trial's condition
        E.ptype = Design(condition, 2);        % based on the condition code, find the current trial's probe type in the Design matrix
        cueing = Cueing(Design(condition, 1));  % based on the condition code, find the current trial's cueing condition in the Design matrix
        output = Model(P, setsize, cueing);     % run the model, get the output structure
        response(ConditionCount(condition), condition) = output.response(1);   % read out the response into the response matrix
        rt(ConditionCount(condition), condition) = output.rt;                  % read out the RT into the RT matrix
        delta = output.response(2);                              % the second entry of response is the degree of change of the probe
        bin = find(delta < binBounds, 1);                        % find the right bin for the degree of change
        Bin(ConditionCount(condition), condition) = bin;
    end
    
    Big = Bin>=8 | Bin==1;
    for condition = 1:12
        ptype = Design(condition, 2);
        ccond = Design(condition, 1);
        Pyes(id, ptype, ccond) = mean(2-response(:,condition));  % Yes-response is coded as 1, No as 2 --> add 1 for every "Yes" and 0 for every "No"
        PyesBig(id, ptype, cueing) = mean(2-response(Big(:,condition),condition));  % For big changes only
        if (ptype == 1)
            PC(id, ptype, ccond) = Pyes(id, ptype, ccond);
            PCbig(id, ptype, cueing) = PyesBig(id, ptype, cueing);
        else
            PC(id, ptype, ccond) = 1-Pyes(id, ptype, ccond);   % for change probes, accuracy = 1-P(yes)
            PCbig(id, ptype, cueing) = 1-PyesBig(id, ptype, cueing);
        end
        RT(id, ptype, ccond) = mean(rt(:,condition));
    end
    
end  % for ID

% Plots

% Proportion correct and RT by cueing condition
PreFigure;
subplot(1,2,1);
plotvector = squeeze(mean(mean(PC, 2), 1));
plot(1:3, plotvector);
PostFigure([0.5, 3.5, 0.5, 1], 'Cue Condition', 'P(correct)');
xticklabels({'No Cue', 'RS Delay', 'Retro-Cue'});
subplot(1,2,2);
plotvector = squeeze(mean(mean(RT, 2), 1));
plot(1:3, plotvector);
PostFigure([0.5, 3.5, 0.5, 2], 'Cue Condition', 'RT(s)');
xticklabels({'No Cue', 'RS Delay', 'Retro-Cue'});

% Proportion correct, separately for hits and FAs
PreFigure;
subplot(1,3,1);
plot(1:3, squeeze(mean(mean(PC(:,1,:), 2), 1)) );
PostFigure([0.5, 3.5, 0.5, 1], 'Cue Condition', 'P(hit)');
xticklabels({'No Cue', 'RS Delay', 'Retro-Cue'});
subplot(1,3,2);
plot(1:3, squeeze(mean(PC(:,2,:), 1)) );
PostFigure([0.5, 3.5, 0.5, 1], 'Cue Condition', 'P(CR new)');
xticklabels({'No Cue', 'RS Delay', 'Retro-Cue'});
subplot(1,3,3);
plot(1:3, squeeze(mean(PC(:,3,:), 1)) );
PostFigure([0.5, 3.5, 0.5, 1], 'Cue Condition', 'P(CR intrus)');
xticklabels({'No Cue', 'RS Delay', 'Retro-Cue'});

%%% Save results
if E.saveResults == 1
    filename = ['IMSim.DelayRS', mat2str(E.material), '.dat'];
    fid = fopen(filename, 'w');
    for id = 1:E.nsubj
        for ptype = 1:3
            for cueing = 1:3
                fprintf(fid, '%d %d %d %d %d %d %d %d \n', id, setsize, ptype, cueing, PC(id, ptype, cueing), RT(id, ptype, cueing), PCbig(id, ptype, cueing), RTbig(id, ptype, cueing));
            end
        end
    end
    fclose(fid);
end
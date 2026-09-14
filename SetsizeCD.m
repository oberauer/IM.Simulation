function [] = SetsizeCD(Model)
% Simulation of Interaction of retro-cueing with set-size in Change Detection
% Experiment of Souza & Oberauer (2014, APP), varying set size (1 to 6) and
% retro-cueing (no cue vs. valid cue). 
% Not-matching probes are sampled from a uniform distribution of the color
% wheel (as in Lin & Oberauer, 2022). 

global P
global C
global E

E.PreRetro = 2;  % this is all retro-cue
E.cuevalidity = 1;
E.maxsetsize = 6; 

C.nstim = 360; 
IMprepareRecog; % set up criterion for expected size of change 

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff; 

% Initializing some container matrices
Pyes = zeros(E.nsubj, E.maxsetsize, 3, 2);  % Probability of saying "Yes" (="Same") for each subject, setsize, probe type, and cueing condition
PyesBig = zeros(E.nsubj, E.maxsetsize, 3, 2);  % Probability of saying "Yes" (="Same") for each subject, setsize, probe type, and cueing condition, for big (or 0) changes only
PC = zeros(E.nsubj, E.maxsetsize, 3, 2);    % Proportion correct
PCbig = zeros(E.nsubj, E.maxsetsize, 3, 2);    % Proportion correct
RT = zeros(E.nsubj, E.maxsetsize, 3, 2);    % Response time
RTbig = zeros(E.nsubj, E.maxsetsize, 3, 2);    % Response time
binBounds = [1, 15:15:180];                 % Bin boundaries for binning the degrees of change among change trials
binBounds(end) = binBounds(end)+0.5; 
nbins = length(binBounds);                  % number of bins
PyesXDelta = zeros(E.nsubj, E.maxsetsize, nbins);  % Probability of saying "Yes" for each bin of degree of change
PyesXDcount = zeros(E.nsubj, E.maxsetsize, nbins);  % Probability of saying "Yes" - counters
LDsame = NaN(E.nsubj, E.maxsetsize);  % Length of the region of similarity between retrieved feature and probe for which "Same" will be responded

for id = 1:E.nsubj
    
    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']); 
    end
    
    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);
    
    for setsize = 1:E.maxsetsize
        
        % Setting up the experimental design
        Cueing = [ones(1,4), 2*ones(1,4)];  % 2 cueing conditions to be crossed with 4 probe-type conditions
        Ptype = [1 1 2 3];  % 2 x positive, 1 x new, 1 x intrusion
        Design = [Cueing', repmat(Ptype', 2, 1)];   % Combining the 2 IV to 2x4=8 design cells
        response = zeros(E.ntrials, 8);             % responses in all trials for the 8 design cells
        rt = zeros(E.ntrials, 8);                   % response times
        Bin = zeros(E.ntrials, 8);                  % bins for changes of probes
        Conditionvector = repmat(1:8, 1, E.ntrials);  % vector of conditions (design cells) for the 8 x ntrials trials
        Conditionvector = Conditionvector(randperm(length(Conditionvector)));  % shuffle the order of design cells
        ConditionCount = zeros(1,8);                                           % running counter of trials in each condition
        LengthDSame = zeros(1,8*E.ntrials);                                    % length of region of similarity for which "same" is responded
        for trial = 1:(8*E.ntrials)
            condition = Conditionvector(trial);         % pick the condition of this trial
            ConditionCount(condition) = ConditionCount(condition) + 1;  % increment trial count for the current trial's condition
            E.ptype = Design(condition, 2);             % determine the probetype from the design matrix
            if (setsize == 1 && Design(condition, 2) == 3), E.ptype = 2; end % for set size 1, there are no intrusion probes
            cueing = Design(condition, 1);              % determine the cueing condition from the design matrix
            output = Model(P, setsize, cueing);   % run model on 1 trial, returns predictions (output is a structure with lots of variables in it)
            response(ConditionCount(condition), condition) = output.response(1,:);  % the first entry of response is the actual response
            rt(ConditionCount(condition), condition) = output.rt;    % response time
            %deltaOld = output.response(3);                              % the third entry of response is the degree of change of the probe
            delta = abs(wrap(C.feature(output.F(1))-C.feature(output.probeIdx), 180)); % size of change of the probe relative to the target feature
            bin = find(delta < binBounds, 1);                        % find the right bin for the degree of change
            PyesXDelta(id, setsize, bin) = PyesXDelta(id, setsize, bin) + (2-output.response(1)); % response is coded Yes=1/No=2, so here we add 1 for Yes, and 0 for No
            PyesXDcount(id, setsize, bin) = PyesXDcount(id, setsize, bin) + 1;         % counting up the number of observations in each bin
            LengthDSame(trial) = sum(output.SameRange);   % SameRange is a vector of 1 in the "same" range and 0 in the "change" range, so their sum is the length of the "same" range
            Bin(ConditionCount(condition), condition) = bin;
        end
        
        % now loop over the 8 design cells to read out the summary
        % statistics of simulated data in each cell
        Big = Bin>=8 | Bin==1;
        for condition = 1:8
            ptype = Design(condition,2);
            cueing = Design(condition,1);
            Pyes(id, setsize, ptype, cueing) = mean(2-response(:,condition));  % Yes/No: response = 1/2
            PyesBig(id, setsize, ptype, cueing) = mean(2-response(Big(:,condition),condition));  % For big changes only
            if (ptype == 1)
                PC(id, setsize, ptype, cueing) = Pyes(id, setsize, ptype, cueing);
                PCbig(id, setsize, ptype, cueing) = PyesBig(id, setsize, ptype, cueing);
            else
                PC(id, setsize, ptype, cueing) = 1-Pyes(id, setsize, ptype, cueing);
                PCbig(id, setsize, ptype, cueing) = 1-PyesBig(id, setsize, ptype, cueing);
            end
            RT(id, setsize, ptype, cueing) = mean(rt(:,condition))./1000;
            RTbig(id, setsize, ptype, cueing) = mean(rt(Big(:,condition),condition))./1000;
        end
        
        LDsame(id, setsize) = mean(LengthDSame);
        
        disp(['      ID      Setsize   PC        RT']);
        disp([id, setsize, mean(mean(PC(id, setsize, :, :))), mean(mean(RT(id, setsize, :, :)))]);
        
    end % for setsize
    
end  % for ID

%%% Plots

% Proportion correct as a function of set size and cueing condition

PreFigure;
subplot(1,2,1);
plotvector(:,1) = squeeze(mean(mean(PC(:,:,:,1), 3), 1));  % no-cue condition
plotvector(:,2) = squeeze(mean(mean(PC(:,:,:,2), 3), 1));  % cue condition
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0.5, 1], 'Set Size', 'P(correct)', 'All changes', {'No Cue', 'Cue'});
subplot(1,2,2);
plotvector(:,1) = squeeze(mean(mean(RT(:,:,:,1), 3), 1));  % no-cue condition
plotvector(:,2) = squeeze(mean(mean(RT(:,:,:,2), 3), 1));  % cue condition
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0, 2], 'Set Size', 'RT', 'All changes', {'No Cue', 'Cue'});

% Accuracy, separately for Hits and FAs

PreFigure;
subplot(1,3,1);
plotvector(:,1) = squeeze(mean(mean(PC(:,:,1,1),3), 1));  % positive probes, no-cue condition
plotvector(:,2) = squeeze(mean(mean(PC(:,:,1,2),3), 1));  % positive probes, cue condition
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0.5, 1], 'Set Size', 'P(hit)', [], {'No Cue', 'Cue'});
subplot(1,3,2);
plotvector(:,1) = squeeze(mean(PC(:,:,2,1), 1));  % new probes, no-cue condition
plotvector(:,2) = squeeze(mean(PC(:,:,2,2), 1));  % new probes, cue condition
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0.5, 1], 'Set Size', 'P(CR New)', [], {'No Cue', 'Cue'});
subplot(1,3,3);
plotvector(:,1) = squeeze(mean(PC(:,:,3,1), 1));  % intrusion probes, no-cue condition
plotvector(:,2) = squeeze(mean(PC(:,:,3,2), 1));  % intrusion probes, cue condition
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0.5, 1], 'Set Size', 'P(CR Intrus)', [], {'No Cue', 'Cue'});

% P(yes) by distance between probe and target feature
PreFigure([], [], 2);
for setsize = 1:E.maxsetsize
    pyes(setsize,:) = squeeze(mean(PyesXDelta(:,setsize,:),1)./mean(PyesXDcount(:,setsize,:),1));
end
plot([0, binBounds(2:end)-7.5], pyes');
PostFigure([0, 180, 0, 1], 'D(probe, target)', 'P(yes)', [], vec2legend(1:setsize));

% For big changes only:
% Proportion correct as a function of set size and cueing condition

PreFigure;
subplot(1,2,1);
plotvector(:,1) = squeeze(mean(mean(PCbig(:,:,:,1), 3), 1));  % no-cue condition
plotvector(:,2) = squeeze(mean(mean(PCbig(:,:,:,2), 3), 1));  % cue condition
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0.5, 1], 'Set Size', 'P(correct)', 'Big changes only', {'No Cue', 'Cue'});
subplot(1,2,2);
plotvector(:,1) = squeeze(mean(mean(RTbig(:,:,:,1), 3), 1));  % no-cue condition
plotvector(:,2) = squeeze(mean(mean(RTbig(:,:,:,2), 3), 1));  % cue condition
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0, 2], 'Set Size', 'RT', 'Big changes only', {'No Cue', 'Cue'});

%%% Save results
if E.saveResults == 1
    filename = ['IMSim.SetsizeCD', mat2str(E.material), '.dat'];
    fid = fopen(filename, 'w');
    for id = 1:E.nsubj
        for setsize = 1:E.maxsetsize
            for ptype = 1:3
                for cueing = 1:2
                    fprintf(fid, '%d %d %d %d %d %d %d %d \n', id, setsize, ptype, cueing, PC(id, setsize, ptype, cueing), RT(id, setsize, ptype, cueing), PCbig(id, setsize, ptype, cueing), RTbig(id, setsize, ptype, cueing));
                end
            end
        end
    end
    fclose(fid);
    filename = ['IMSim.SetsizeCD.PYesByBin.', mat2str(E.material), '.dat'];
    fid = fopen(filename, 'w');
    for setsize = 1:E.maxsetsize
        for bin = 1:nbins, fprintf(fid, '%d ', pyes(setsize, bin)); end
        fprintf(fid, '\n');
    end
    fclose(fid);
end





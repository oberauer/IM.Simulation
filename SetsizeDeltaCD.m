function D = SetsizeDeltaCD(Model)

% Not-matching probes are sampled from a uniform distribution of the color
% wheel (as in Lin & Oberauer, 2022). 

global P
global C
global E

E.test = 2;      % change detection
E.PreRetro = 2;  % this is all retro-cue
E.maxsetsize = 6; 

likSame = VonMisesN(C.x, pi, P.kappacrit);  % estimated likelihood of "same" trials using the meta-cognitive estimate of feature precision
likChange = 1./360;  % estimated likelihood of "change" trials: all 360 values equally likely as change probes
E.SameRange = likSame > likChange;  % range of retrieved feature values that are similar enough to the probe to say "same"


% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff; 

% Set up the design
cueing = 1;
Ptype = [1, 1, 2, 3];   % levels of the probetype variable 2 x positive, 1 x new, 1 x intrusion
Design = fullfact(4);   % probetype x cueing
nCells = size(Design, 1);      % number of design cells (crossing cueing with probetype)

% Initializing some container matrices
Pyes = zeros(E.nsubj, E.maxsetsize, 3);  % Probability of saying "Yes" (="Same") for each subject, setsize, probe type
PC = zeros(E.nsubj, E.maxsetsize, 3);    % Proportion correct
RT = zeros(E.nsubj, E.maxsetsize, 3);    % Response time
binBounds = [1, 15:15:180];                 % Bin boundaries for binning the degrees of change among change trials
binBounds(end) = binBounds(end)+0.5; 
nbins = length(binBounds);                  % number of bins
PyesXDelta = zeros(E.nsubj, E.maxsetsize, nbins);  % Probability of saying "Yes" for each bin of degree of change
PyesXDcount = zeros(E.nsubj, E.maxsetsize, nbins);  % Probability of saying "Yes" - counters
PyesXDelta_Alt = zeros(E.nsubj, E.maxsetsize, nbins); 

PyesXDelta_New = zeros(E.nsubj, E.maxsetsize, nbins); 
PyesXDelta_Intrus = zeros(E.nsubj, E.maxsetsize, nbins); 
PyesXDelta_NewCount = zeros(E.nsubj, E.maxsetsize, nbins); 
PyesXDelta_IntrusCount = zeros(E.nsubj, E.maxsetsize, nbins); 

LDsame = NaN(E.nsubj, E.maxsetsize);  % Length of the region of similarity between retrieved feature and probe for which "Same" will be responded

for id = 1:E.nsubj
    
    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']); 
    end
    
    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);

    fdistance = zeros(E.ntrials, E.maxsetsize);  % distance (target, response)
    
    for setsize = 1:E.maxsetsize
        
        response = zeros(E.ntrials, nCells);             % responses in all trials for the 4 design cells
        rt = zeros(E.ntrials, nCells);                   % response times
        Bin = zeros(E.ntrials, nCells);                  % bins for changes of probes
        Conditionvector = repmat(1:nCells, 1, E.ntrials);  % vector of conditions (design cells) for the nCells x ntrials trials
        Conditionvector = Conditionvector(randperm(length(Conditionvector)));  % shuffle the order of design cells
        ConditionCount = zeros(1,nCells);                                           % running counter of trials in each condition
        LengthDSame = zeros(1,nCells*E.ntrials);                                    % length of region of similarity for which "same" is responded
        for trial = 1:(nCells*E.ntrials)
            condition = Conditionvector(trial);         % pick the condition of this trial
            ConditionCount(condition) = ConditionCount(condition) + 1;  % increment trial count for the current trial's condition
            E.ptype = Ptype(Design(condition, 1));             % determine the probetype from the design matrix
            if (setsize == 1 && Design(condition, 1) == 3), E.ptype = 2; end % for set size 1, there are no intrusion probes
            output = Model(P, setsize, cueing);   % run model on 1 trial, returns predictions (output is a structure with lots of variables in it)
            response(ConditionCount(condition), condition) = output.response(1,:);  % the first entry of response is the actual response
            rt(ConditionCount(condition), condition) = output.rt;    % response time
            delta = abs(wrap(C.feature(output.F(1))-C.feature(output.probeIdx), 180)); % size of change of the probe relative to the target feature
            bin = find(delta < binBounds, 1);                        % find the right bin for the degree of change
            PyesXDelta(id, setsize, bin) = PyesXDelta(id, setsize, bin) + (2-output.response(1)); % response is coded Yes=1/No=2, so here we add 1 for Yes, and 0 for No
            PyesXDcount(id, setsize, bin) = PyesXDcount(id, setsize, bin) + 1;         % counting up the number of observations in each bin
            LengthDSame(trial) = sum(output.SameRange);   % SameRange is a vector of 1 in the "same" range and 0 in the "change" range, so their sum is the length of the "same" range
            Bin(ConditionCount(condition), condition) = bin;
            fdistance(trial, setsize) = wrap(output.response(2)-C.feature(output.probeIdx), 180);  % response(2) is the retrieved feature: how far is it from the probe?
            PyesXDelta_Alt(id, setsize, bin) = PyesXDelta_Alt(id, setsize, bin) + double(abs(fdistance(trial, setsize)) < LengthDSame(trial)./2); % if retrieved feature lies within "Same" interval, "yes" response is expected
            if E.ptype == 1 || E.ptype == 2
                PyesXDelta_New(id, setsize, bin) = PyesXDelta_New(id, setsize, bin) + (2-output.response(1)); % response is coded Yes=1/No=2, so here we add 1 for Yes, and 0 for No
                PyesXDelta_NewCount(id, setsize, bin) = PyesXDelta_NewCount(id, setsize, bin) + 1;         % counting up the number of observations in each bin
            end
            if E.ptype == 1 || E.ptype == 3
                PyesXDelta_Intrus(id, setsize, bin) = PyesXDelta_Intrus(id, setsize, bin) + (2-output.response(1)); % response is coded Yes=1/No=2, so here we add 1 for Yes, and 0 for No
                PyesXDelta_IntrusCount(id, setsize, bin) = PyesXDelta_IntrusCount(id, setsize, bin) + 1;         % counting up the number of observations in each bin
            end
        end
        
        % now loop over the 4 design cells to read out the summary statistics of simulated data in each cell
        for condition = 1:4
            ptype = Ptype(Design(condition,2));
            cueing = Design(condition,1);
            Pyes(id, setsize, ptype) = Pyes(id, setsize, ptype) + mean(2-response(:,condition));  % Yes/No: response = 1/2
            if (ptype == 1)
                PC(id, setsize, ptype) = PC(id, setsize, ptype) + Pyes(id, setsize, ptype)./2;  % divide by 2 because there are 2 conditions for positive probes
            else
                PC(id, setsize, ptype) = 1-Pyes(id, setsize, ptype);
            end
            RT(id, setsize, ptype) = mean(rt(:,condition));
        end
        
        LDsame(id, setsize) = mean(LengthDSame);
        
        disp(['      ID      Setsize   PC        RT']);
        disp([id, setsize, mean(mean(PC(id, setsize, :))), mean(mean(RT(id, setsize, :)))]);
        
    end % for setsize
    
end  % for ID

%%% Plots

% Proportion correct as a function of set size and cueing condition

PreFigure;
subplot(1,2,1);
plotvector = squeeze(mean(PC, 1));  % average over subjects
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0.5, 1], 'Set Size', 'P(correct)', [], {'Positive', 'New', 'Intrusion'});
subplot(1,2,2);
plotvector = squeeze(mean(RT, 1));  % average over subjects
plot(1:E.maxsetsize, plotvector);
PostFigure([0.5, E.maxsetsize+0.5, 0, 2], 'Set Size', 'RT', [], {'Positive', 'New', 'Intrusion'});

% P(yes) by distance between probe and target feature
PreFigure([], [], 2);
for setsize = 1:E.maxsetsize
    pyes(setsize,:) = squeeze(mean(PyesXDelta(:,setsize,:),1)./mean(PyesXDcount(:,setsize,:),1));
end
plot([0, binBounds(2:end)-7.5], pyes');
PostFigure([0, 180, 0, 1], 'D(probe, target)', 'P(yes)', [], vec2legend(1:setsize));

% P(yes) by distance between probe and target feature, separately for new and intrusion lures
PreFigure([], [], 2);
subplot(1,2,1);
for setsize = 1:E.maxsetsize
    pyes(setsize,:) = squeeze(mean(PyesXDelta_New(:,setsize,:),1)./mean(PyesXDelta_NewCount(:,setsize,:),1));
end
plot([0, binBounds(2:end)-7.5], pyes');
PostFigure([0, 180, 0, 1], 'D(probe, target)', 'P(yes)', 'New', vec2legend(1:setsize));
subplot(1,2,2);
for setsize = 1:E.maxsetsize
    pyes(setsize,:) = squeeze(mean(PyesXDelta_Intrus(:,setsize,:),1)./mean(PyesXDelta_IntrusCount(:,setsize,:),1));
end
plot([0, binBounds(2:end)-7.5], pyes');
PostFigure([0, 180, 0, 1], 'D(probe, target)', 'P(yes)', 'Intrus', vec2legend(1:setsize));

D.PC = PC;
D.Pyes = Pyes;
D.RT = RT;

%%% Save results
if E.saveResults == 1
    filename = ['IMSim.SetsizeCD', mat2str(E.material), '.dat'];
    fid = fopen(filename, 'w');
    for id = 1:E.nsubj
        for setsize = 1:E.maxsetsize
            for ptype = 1:3
                    fprintf(fid, '%d %d %d %d %d \n', id, setsize, ptype, PC(id, setsize, ptype), RT(id, setsize, ptype));
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





function [] = MultiCueABA(Model)
% Simulation of retro-cue effect in ABA vs. CBA sequence (Rerko et al., 2013, Exp. 3)

global P
global E
global C

E.PreRetro = 2;  % this is all retro-cue
E.material = 2;
setsize = 6;
E.cuevalidity = 1/setsize;  % analogous to the refreshing experiments, because the 1st and 2nd cue are often not valid

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

IMprepareRecog;

% generate parameters with individual differences
ParX = CreateIndDiff;

% Initializing some container matrices
Pyes = zeros(E.nsubj, 3, 3);  % Probability of saying "Yes" (="Same") for each subject, probe type, and cueing condition
PC = zeros(E.nsubj, 3, 3);    % Proportion correct
RT = zeros(E.nsubj, 3, 3);    % Response time
Error = zeros(E.nsubj, 3, 3); % Continuous-reproduction error

for test = 1:2
    E.test = test;

    for id = 1:E.nsubj

        % extract parameter values for each subject - for those parameters that vary between subjects
        for ii = 1:length(C.indVar)
            eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
        end

        % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
        CreateStimuli;
        CreateMapping(E.calibrateAmp==2);

        % Setting up the experimental design
        Cueing = [ones(1,4), 2*ones(1,4), 3*ones(1,4)];  % 3 cueing conditions to be crossed with 4 probe-type conditions
        Ptype = [1 1 2 3];  % 2 x positive, 1 x new, 1 x intrusion
        Design = [Cueing', repmat(Ptype', 3, 1)];   % Combining the 2 IV to 2x4=8 design cells
        response = zeros(E.ntrials, 12);             % responses in all trials for the 12 design cells
        rt = zeros(E.ntrials, 12);                   % response times
        error = zeros(E.ntrials, 12);                % continuous-reproduction errors
        Conditionvector = repmat(1:12, 1, E.ntrials);  % vector of conditions (design cells) for the 12 x ntrials trials
        Conditionvector = Conditionvector(randperm(length(Conditionvector)));  % shuffle the order of design cells
        ConditionCount = zeros(1,12);                                           % running counter of trials in each condition
        
        for trial = 1:(12*E.ntrials)

            condition = Conditionvector(trial);         % pick the condition of this trial
            ConditionCount(condition) = ConditionCount(condition) + 1;  % increment trial count for the current trial's condition
            E.ptype = Design(condition, 2);             % determine the probetype from the design matrix
            if (setsize == 1 && Design(condition, 2) == 3), E.ptype = 2; end % for set size 1, there are no intrusion probes
            cueingcond = Design(condition, 1);              % determine the cueing condition from the design matrix
            cueing = 1 + 3*(cueingcond>1); % cueingcond 1 = no cue, 2 & 3 = refreshing --> cueing = 4
            if cueingcond == 2, C.RefSequence(1).seq = [1+randperm(setsize-1, 2), 1]; end   % CBA
            if cueingcond == 3, C.RefSequence(1).seq = [1, 1+randperm(setsize-1, 1), 1]; end  % ABA

            if E.test == 2
                output = Model(P, setsize, cueing);   % run model on 1 trial, returns predictions (output is a structure with lots of variables in it)
                response(ConditionCount(condition), condition) = output.response(1,:);  % the first entry of response is the actual response
                %disp([cueing, E.ptype, output.response(1,:)]);
                rt(ConditionCount(condition), condition) = output.rt;    % response time
            end

            if E.test == 1
                % continuous-reproduction test
                outputR = Model(P, setsize, cueing);
                featurestep = floor(C.nc/C.nstim);
                targetDeg = outputR.F(1)*featurestep; % translate the target feature F (1 to 12) into the angle in degrees
                fdistance = wrap(outputR.response-targetDeg, 180);   %calculate distance between response and true feature in feature space (degrees!)
                error(ConditionCount(condition), condition) = abs(fdistance);
                E.test = 2;
            end

        end

        if E.test == 2
            % now loop over the 12 design cells to read out the summary
            % statistics of simulated data in each cell
            for condition = 1:12
                ptype = Design(condition,2);
                cueingcond = Design(condition,1);
                Pyes(id, ptype, cueingcond) = mean(2-response(:,condition));  % Yes/No: response = 1/2
                if (ptype == 1)
                    PC(id, ptype, cueingcond) = PC(id, ptype, cueingcond) + 0.5*Pyes(id, ptype, cueingcond);  % there are 2 conditions per cell with positive probes
                else
                    PC(id, ptype, cueingcond) = 1-Pyes(id, ptype, cueingcond);
                end
                RT(id, ptype, cueingcond) = mean(rt(:,condition));
                Error(id, ptype, cueingcond) = mean(error(:,condition));
            end
        end

        disp(['ID = ', mat2str(id)]);

    end  % for ID

    %%% Plots

    % Proportion correct and RT as a function of cueing condition

    if E.test == 2
        Legend = {'Positive', 'New', 'Intrusion'};
        PreFigure;
        subplot(1,2,1);
        plot(1:3, squeeze(mean(PC, 1))');
        PostFigure([0.5, 3.5, 0, 1], 'Cueing Condition', 'P(correct)', [], Legend);
        xticklabels({'None','CBA','ABA'});
        subplot(1,2,2);
        plot(1:3, squeeze(mean(RT, 1))');
        PostFigure([0.5, 3.5, 0, 2], 'Cueing Condition', 'RT', [], Legend);
        xticklabels({'None','CBA','ABA'});
    end

    if E.test == 1
        PreFigure;
        plot(1:3, squeeze(mean(Error, 1))');
        PostFigure([0.5, 3.5, 0, 90], 'Cueing Condition', 'Reproduction Error', [], Legend);
        xticks(1:3);
        xticklabels({'None','CBA','ABA'});
    end

end

%%% Save results
if E.saveResults == 1
    filename = ['IMSim.MultiCueABA.dat'];
    fid = fopen(filename, 'w');
    for id = 1:E.nsubj
        for ptype = 1:3
            for cueing = 1:3
                fprintf(fid, '%d %d %d %d %d \n', id, ptype, cueing, PC(id, ptype, cueing), RT(id, ptype, cueing));
            end
        end
    end
    fclose(fid);
end





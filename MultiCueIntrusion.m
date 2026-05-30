function [] = MultiCueIntrusion(Model)
% Simulation of two successive cues in Change Detection, with intrusion
% probes matching the first-cued item (Rerko & Oberauer, 2013, Exp. 2)

global P
global C
global E

E.PreRetro = 2;  % this is all retro-cue
E.material = 2; 
E.cuesequence = [0,1];      % successive cue targets, with "0" for "any non-target chosen at random", "1" for target, and any number > 1 for the specific non-target indexed by that number
E.cuevalidity = 0.5;        % In multi-cue paradigm, predictability of number of cues appears to have no effect, so this must be set low for all cues
IMprepareRecog;

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff; 

% Initializing some container matrices
Pyes = zeros(E.nsubj, 5);  % Probability of saying "Yes" (="Same") for each subject and probe type (pos, new, intrustion neighbor, intrusion non-neighbor, intrusion cued)
PC = zeros(E.nsubj, 5);    % Proportion correct
RT = zeros(E.nsubj, 5);    % Response time

for id = 1:E.nsubj
    
        % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']); 
    end
    
    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);
    
    ConditionCount = ones(1,5);
    response = NaN(E.ntrials,5);
    rt = NaN(E.ntrials,5);
    
    for ptype = 1:3
        
        E.ptype = ptype;
        
        for trial = 1:E.ntrials
            output = Model(P, 6, 6);   % run model on 1 trial, returns predictions (output is a structure with lots of variables in it); setsize = 6, cueing = 6
            probetype = ptype;
            if ptype == 3
                if abs(output.L(output.probeIdx) - output.L(1)) == 1, probetype = 3; else, probetype = 4; end % neighbor vs. non-neighbor intrusion
                if output.probeIdx == output.cueIdx(1), probetype = 5; end % first-cued item becomes intrusion probe
            end
            response(ConditionCount(probetype), probetype) = output.response(1);  % the first entry of response is the actual response
            rt(ConditionCount(probetype), probetype) = output.rt;    % response time
            ConditionCount(probetype) = ConditionCount(probetype) + 1;
        end
        
    end % for ptype
    
    for probetype = 1:5
        Pyes(id, probetype) = nanmean(2-response(:,probetype));  % Yes/No: response = 1/2
        if (probetype == 1)
            PC(id, probetype) = Pyes(id, probetype);
        else
            PC(id, probetype) = 1-Pyes(id, probetype);
        end
        RT(id, probetype) = nanmean(rt(:,probetype));
    end
    
end  % for ID

%%% Plots

% Proportion correct and RT as function of probe type

xTicks = {'Pos', 'New', 'Intrus-Far', 'Intrus-Near', 'Intrus-Cued'};
PreFigure;
subplot(1,2,1);
plot(1:5, mean(PC));
PostFigure([0.5, 5.5, 0.5, 1], 'Probe Type', 'P(correct)');
subplot(1,2,2);
plot(1:5, mean(RT));
PostFigure([0.5, 5.5, 0, 1.1*max(mean(RT))], 'Probe Type', 'RT');
xticklabels(xTicks);

%%% Save results
if E.saveResults == 1
    fid = fopen(['IMSim.MultiCueIntrusion.dat'], 'w');
    for id = 1:E.nsubj
        for ptype = 1:5
            fprintf(fid, '%d %d %d %d \n', id, ptype, PC(id, ptype), RT(id, ptype));
        end
    end
    fclose(fid);
end


function D = Reloading(Model)
% Simulation of Reloading effects: 
% Experiment 2 in Souza, Rerko, & Oberauer (2014, JEP:HPP): Participants encode a
% first array (set size 2 vs. 4), followed by a retro cue (or no cue). Then
% they encode a second array (set size 2 vs. 4), and are then tested on the
% first, then the second, array. 

global P
global E
global C

E.test = 2;      % change detection
E.PreRetro = 2;  % this is all retro-cue
E.cuevalidity = 1;
E.material = 2; 
E.mask = 3;
E.ntrials = round(E.ntrials/12);  % because we're running 32 design cells
IMprepareRecog;

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff; 

% setting up containers for results
Pyes1 = zeros(E.nsubj, 2, 3, 2);  % Prob. of saying "Yes" (="Same") for the item tested in the first array
Pyes2 = zeros(E.nsubj, 2, 3, 2);  % P(yes) for second array, by subject, set size, probe type, and cueing condition
PC1 = zeros(E.nsubj, 2, 3, 2);    % Proportion correct
PC2 = zeros(E.nsubj, 2, 3, 2);
RT1 = zeros(E.nsubj, 2, 3, 2);    % Response time
RT2 = zeros(E.nsubj, 2, 3, 2);
GRecord = zeros(2,2); 
GRecordCount = zeros(2,2); 
BStrength = zeros(2,2);

Cueing = 1:2;           % levels of the cueing variable
Ptype = [1, 1, 2, 3];   % levels of the probetype variable 2 x positive, 1 x new, 1 x intrusion
Design = fullfact([2, 4, 4]);  % cueing, probetype (1st array), probetype (2nd array)
nCells = size(Design, 1);      % number of design cells (crossing cueing, probetype array 1, and probetype array 2)

for id = 1:E.nsubj
    
    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']); 
    end
    
    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli(1);
    CreateMapping(E.calibrateAmp==2);
       
        for ss = 1:2   
            setsize = ss*2;   % set size of array 1 (set size of array 2 is always 4)
            response1 = zeros(E.ntrials, nCells);
            response2 = zeros(E.ntrials, nCells);
            rt1 = zeros(E.ntrials, nCells);
            rt2 = zeros(E.ntrials, nCells);
            Conditionvector = repmat(1:nCells, E.ntrials);
            Conditionvector = Conditionvector(randperm(length(Conditionvector)));  %shuffle condition vector
            ConditionCount = zeros(1,nCells);
            
            for trial = 1:(nCells*E.ntrials)
                
                condition = Conditionvector(trial);
                ConditionCount(condition) = ConditionCount(condition) + 1;  % increment trial count for the current trial's probe type
                E.ptype = Ptype(Design(condition, 2));
                cueing = Cueing(Design(condition, 1));
                E.layout = 3; % all in the first semi-circle
                
                output = Model(P, setsize, cueing);   % first set
                
                response1(ConditionCount(condition), condition) = output.response(1,:);  % get response to the test of the first set (before encoding the second set, so not quite like in the experiment)
                rt1(ConditionCount(condition), condition) = output.rt;           
                Mstate = {output.map, output.wx, output.g};        
                GRecord(ss, cueing) = GRecord(ss, cueing) + mean(output.g); 
                GRecordCount(ss, cueing) = GRecordCount(ss, cueing) + 1; 

                E.ptype = Ptype(Design(condition, 3));
                E.layout = 4; % all in the second semi-circle

                output = Model(P, 4, 1, 0, Mstate);   % second set: setsize=4, neutral cueing, but takes over memory state from first set
                
                response2(ConditionCount(condition), condition) = output.response(1,:);
                rt2(ConditionCount(condition), condition) = output.rt;

                BStrength(ss, cueing) = BStrength(ss, cueing) + mean(output.Bstrength); 

            end
            
            for condition = 1:nCells
                ptype1 = Ptype(Design(condition,2));
                ptype2 = Ptype(Design(condition,3));
                cueing = Cueing(Design(condition,1));
                Pyes1(id, ss, ptype1, cueing) = mean(2-response1(:,condition));  % Yes/No: response = 1/2
                Pyes2(id, ss, ptype2, cueing) = mean(2-response2(:,condition));
                if (ptype1 == 1)
                    PC1(id, ss, ptype1, cueing) = Pyes1(id, ss, ptype1, cueing);
                else
                    PC1(id, ss, ptype1, cueing) = 1-Pyes1(id, ss, ptype1, cueing);
                end
                if (ptype2 == 1)
                    PC2(id, ss, ptype2, cueing) = Pyes2(id, ss, ptype2, cueing);
                else
                    PC2(id, ss, ptype2, cueing) = 1-Pyes2(id, ss, ptype2, cueing);
                end
                RT1(id, ss, ptype1, cueing) = mean(rt1(:,condition));
                RT2(id, ss, ptype2, cueing) = mean(rt2(:,condition));
            end
                        
        disp('     ID  Setsize');
        disp([id, setsize]);

        end % for setsize
           
end  % for ID


meanGRecord = GRecord ./ GRecordCount; 
disp(meanGRecord); 
meanBStrength = BStrength ./ GRecordCount; 
disp(meanBStrength); 

% Plots

PreFigure;
subplot(2,2,1);
plotvector(:,1) = squeeze(mean(mean(PC1(:,:,:,1), 3), 1));  % no-cue condition
plotvector(:,2) = squeeze(mean(mean(PC1(:,:,:,2), 3), 1));  % cue condition
plot([2,4], plotvector);
PostFigure([1, 4.5, 0.5, 1], 'Set Size', 'P(correct)', 'Response 1', {'No Cue', 'Cue'});
subplot(2,2,2);
plotvector(:,1) = squeeze(mean(mean(PC2(:,:,:,1), 3), 1));  % no-cue condition
plotvector(:,2) = squeeze(mean(mean(PC2(:,:,:,2), 3), 1));  % cue condition
plot([2,4], plotvector);
PostFigure([1, 4.5, 0.5, 1], 'Set Size', 'P(correct)', 'Response 2', {'No Cue', 'Cue'});
subplot(2,2,3);
plotvector(:,1) = squeeze(mean(mean(RT1(:,:,:,1), 3), 1));  % no-cue condition
plotvector(:,2) = squeeze(mean(mean(RT1(:,:,:,2), 3), 1));  % cue condition
plot([2,4], plotvector);
PostFigure([1, 4.5, 0, 2], 'Set Size', 'RT (s)', 'Response 1', {'No Cue', 'Cue'});
subplot(2,2,4);
plotvector(:,1) = squeeze(mean(mean(RT2(:,:,:,1), 3), 1));  % no-cue condition
plotvector(:,2) = squeeze(mean(mean(RT2(:,:,:,2), 3), 1));  % cue condition
plot([2,4], plotvector);
PostFigure([1, 4.5, 0, 2], 'Set Size', 'RT (s)', 'Response 2', {'No Cue', 'Cue'});

% Accuracy, separately for Hits and FAs

PreFigure;
subplot(2,3,1);
plotvector(:,1) = squeeze(mean(mean(PC1(:,:,1,1),3), 1));  % positive probes, no-cue condition
plotvector(:,2) = squeeze(mean(mean(PC1(:,:,1,2),3), 1));  % positive probes, cue condition
plot([2,4], plotvector);
PostFigure([1, 4.5, 0.5, 1], 'Set Size', 'P(hit)', 'Response 1', {'No Cue', 'Cue'});
subplot(2,3,2);
plotvector(:,1) = squeeze(mean(PC1(:,:,2,1), 1));  % new probes, no-cue condition
plotvector(:,2) = squeeze(mean(PC1(:,:,2,2), 1));  % new probes, cue condition
plot([2,4], plotvector);
PostFigure([1, 4.5, 0.5, 1], 'Set Size', 'P(CR New)', 'Response 1', {'No Cue', 'Cue'});
subplot(2,3,3);
plotvector(:,1) = squeeze(mean(PC1(:,:,3,1), 1));  % intrusion probes, no-cue condition
plotvector(:,2) = squeeze(mean(PC1(:,:,3,2), 1));  % intrusion probes, cue condition
plot([2,4], plotvector);
PostFigure([1, 4.5, 0.5, 1], 'Set Size', 'P(CR Intrus)', 'Response 1', {'No Cue', 'Cue'});

subplot(2,3,4);
plotvector(:,1) = squeeze(mean(mean(PC2(:,:,1,1),3), 1));  % positive probes, no-cue condition
plotvector(:,2) = squeeze(mean(mean(PC2(:,:,1,2),3), 1));  % positive probes, cue condition
plot([2,4], plotvector);
PostFigure([1, 4.5, 0.5, 1], 'Set Size', 'P(hit)', 'Response 2', {'No Cue', 'Cue'});
subplot(2,3,5);
plotvector(:,1) = squeeze(mean(PC2(:,:,2,1), 1));  % new probes, no-cue condition
plotvector(:,2) = squeeze(mean(PC2(:,:,2,2), 1));  % new probes, cue condition
plot([2,4], plotvector);
PostFigure([1, 4.5, 0.5, 1], 'Set Size', 'P(CR New)', 'Response 2', {'No Cue', 'Cue'});
subplot(2,3,6);
plotvector(:,1) = squeeze(mean(PC2(:,:,3,1), 1));  % intrusion probes, no-cue condition
plotvector(:,2) = squeeze(mean(PC2(:,:,3,2), 1));  % intrusion probes, cue condition
plot([2,4], plotvector);
PostFigure([1, 4.5, 0.5, 1], 'Set Size', 'P(CR Intrus)', 'Response 2', {'No Cue', 'Cue'});

D.PC1 = PC1;
D.PC2 = PC2;
D.RT1 = RT1;
D.RT2 = RT2;


%%% Save results
if E.saveResults == 1
    filename = ['IMSim.Reloading', mat2str(E.material), '.dat'];
    fid = fopen(filename, 'w');
    for id = 1:E.nsubj
        for setsize = 1:2
            for ptype = 1:3
                for cueing = 1:2
                    fprintf(fid, '%d %d %d %d %d %d %d %d \n', id, setsize, ptype, cueing, PC1(id, setsize, ptype, cueing), RT1(id, setsize, ptype, cueing), PC2(id, setsize, ptype, cueing), RT2(id, setsize, ptype, cueing));
                end
            end
        end
    end
    fclose(fid);
end




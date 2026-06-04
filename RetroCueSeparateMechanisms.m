function [] = RetroCueSeparateMechanisms(Model, Mechanisms, Tasks, Cueconds, fitMM)
% Simulation of retro-cue effect with each mechanism (strengthening, removal, visual interference) switched on individually,
% varying task (continuous reproduction vs. change detection)

global P
global E
global C

E.PreRetro = 2;
E.cuevalidity = 1;
setsize = 6;
Titletext = {'Strength. Only', 'Consolid. Only', 'Removal Only', 'Visual Interf. Only', 'FX Only ', 'Headstart Only', 'None'};

IMprepareRecog; % set up criterion for expected size of change 
Ptype = [1 1 2 3];  % 2 x positive, 1 x new, 1 x intrusion
option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

meanAcc = zeros(2, 4, 3);  % task, mechanism, cueing condition

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% keep default values
%cueStrength = P.cueingStrength;
%retroCueConsolid = C.retroCueConsolid;
filter = P.filter; 
removalThreshold = P.removalThreshold;
eraseFX = P.eraseFX;


% generate parameters with individual differences
ParX = CreateIndDiff; 
CueRate = P.cuerate;   % remember for re-instatement (because it is not an individual-differences variable that will be re-instated from ParX)

for task = Tasks
    
    E.test = task;
    E.material = task; % for change detection, use small set of highly distinct stimuli
    
    for mechanism = Mechanisms

        % eliminate all retro-cue mechanisms
        P.cueingStrength = 0; 
        P.removalThreshold = 0; 
        P.filter = zeros(1,3); 
        C.retroCueConsolid = 0; 
        E.CTI(2) = 0; 
        P.eraseFX = 0;
        
        Mdevobs = NaN(E.nsubj, 3);  % mean observed deviation of responses from target feature for subjects, cueing conditions,
        Mrt = NaN(E.nsubj, 3);      % mean RT
        MMSD = NaN(E.nsubj, 3);     % Mixture Model SD parameter
        MMguessing = NaN(E.nsubj, 3);  % Mixture Model Guessing parameter
        MMtranspos = NaN(E.nsubj, 3);  % Mixture Model Transposition parameter (swap error proportion)
        MMcwattraction = NaN(E.nsubj, 3);  % Mixture Model colorwheel-attraction strength parameter
        Pcorrect = NaN(E.nsubj, 3, 3);       % mean accuracy in change detection for 3 cueing conditions and 3 probe types
        
        % for each set-size level, for each trial, generate a vector of 360 color
        % values coding the colors on the wheel
        [aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x E.maxsetsize x [1:360]
        
        for id = 1:E.nsubj
            
            % extract parameter values for each subject - for those parameters that vary between subjects
            for ii = 1:length(C.indVar)
                eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
            end

            % re-introduce individual mechanisms
            if mechanism == 1, P.cueingStrength = 1; end  % leave only strengthening
            if mechanism == 2, C.retroCueConsolid = 1; end  % leave ony consolidation for retrieval
            if mechanism == 3, P.removalThreshold = removalThreshold; end  % leave only removal
            if mechanism == 4, P.filter = filter; end  % leave only visual interference
            if mechanism == 5, E.CTI(2) = 1; P.eraseFX = eraseFX; end  % leave only head start for read-out from FX
            if mechanism == 6, E.CTI(2) = 1; end  % leave only headstart
            
            % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
            CreateStimuli;
            CreateMapping(E.calibrateAmp==2);
            
            % Initialize container vectors
            fdistance = zeros(1,E.ntrials);  % feature distance between response and target
            rt = zeros(1,E.ntrials);         % response time
            Probedpos = zeros(E.ntrials,1);  % Number of the tested (probed) spatial position
            Pangle = zeros(E.ntrials,setsize);  % spatial angles of item positions in the array
            Cangle = zeros(E.ntrials,setsize+1); % color angles in the color wheel
            Targ = zeros(E.ntrials,1);                % Target
            Resp = zeros(E.ntrials,1);                % Response
            Setsize = zeros(E.ntrials,1);
            Probetype = zeros(E.ntrials,1);             % Probe type for CD
            
            for cueing = Cueconds  % neutral, valid, invalid
                
                for trial = 1:E.ntrials
                    
                    if task == 2
                        E.ptype = Ptype(mod(trial,4)+1);
                    end
                    
                    output = Model(P, setsize, cueing);   % here the model is run!
                    
%                     Array(tcount,:) = output.F(1,1:setsize);  % record of the array on this trial
%                     Target(tcount) = output.F(1,1);    % target feature
                    
                    if task == 1
                        %Response(tcount) = output.response;   % response feature
                        fdistance(trial) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
                    end
                    if task == 2
                        Resp(trial) = output.response(1,:);  % the first entry of response is the actual response
                        Probetype(trial) = E.ptype;
                        %delta = output.response(2);
                    end
                    rt(trial) = output.rt;            % response time
                    
                    %tcount = tcount+1;   % trial counter is incremented
                    
                    if task == 1
                        %collect data for further modeling with Mixture Model
                        Probedpos(trial) = output.L(1);    % probed position
                        Pangle(trial,:) = output.L(1:setsize);  % spatial position angles
                        Cangle(trial,1:(setsize+1)) = [output.F(1,1:setsize), output.CWcolor];  % array feature angles, followed by color-wheel feature angles
                        Targ(trial) = output.F(1,1);           % target feature
                        Resp(trial) = output.response;       % response feature
                        Setsize(trial) = setsize;
                    end
                    
                end
                
                if task == 1
                    Mdevobs(id, cueing) = mean(abs(fdistance));  %mean deviation (averaged over trials)
                end
                if task == 2
                    Pyes = 2-Resp;  % Yes/No: response = 1/2
                    Pcorrect(id, cueing, 1) = mean(Pyes(Probetype==1));
                    for ptype = 2:3
                        Pcorrect(id, cueing, ptype) = mean(1-Pyes(Probetype==ptype));
                    end
                end
                Mrt(id, cueing) = mean(rt);
                
                if fitMM == 1 && task == 1
                    ssData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting
                    
                    % fit Mixture Model
                    startparms = [15, .1, .1, .1];
                    lb = [eps, 0, 0, 0]; ub = [90, 1, 1, 1];
                    npar = 4;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture, 4 -> Souza & Oberauer mixture (iC.ncluding color-wheel attraction)
                    MMloglik = 500000;
                    itercount = 0;
                    while MMloglik > 400000
                        [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, ssData, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                        itercount = itercount + 1;
                    end
                    disp('    id       cueing    mechanism  Deviance  SD(deg)   P(guess)  P(trans)   P(CW)');
                    disp([id, cueing, mechanism, MMloglik/1000, MMparms]);
                    MMSD(id, cueing) = MMparms(1);
                    MMguessing(id, cueing) = MMparms(2);
                    if npar > 2, MMtranspos(id, cueing) = MMparms(3); end
                    if npar > 3, MMcwattraction(id, cueing) = MMparms(4); end
                else
                    disp('    id    task  cueing mechanism ');
                    disp([id, task, cueing, mechanism]);
                end
                
            end % for cueing
            
        end  % for ID
     
        if task == 1, Accuracy = Mdevobs; end
        if task == 2, Accuracy = squeeze(mean(Pcorrect,3)); end  % average over probe types
        
        meanAcc(task, mechanism, :) = mean(Accuracy, 1);
        
        % Plot Mixture Model Parameters over Setsize
        if fitMM && task == 1
            MMPm = 1 - MMtranspos - MMguessing - MMcwattraction;
            PreFigure
            subplot(3,2,1);
            plot(mean(MMSD,1));
            PostFigure([0.5, 3.5, 0, max(mean(MMSD,1))+0.5], 'Cueing Condition', 'Mean SD', Titletext{mechanism});
            subplot(3,2,2);
            plot(mean(MMPm,1));
            PostFigure([0.5, 3.5, 0, 1], 'Cueing Condition', 'Mean P(m)', Titletext{mechanism});
            subplot(3,2,3);
            plot(mean(MMguessing,1));
            PostFigure([0.5, 3.5, 0, 0.5], 'Cueing Condition', 'Mean P(guess)', Titletext{mechanism});
            subplot(3,2,4);
            plot(mean(MMtranspos,1));
            PostFigure([0.5, 3.5, 0, 1], 'Cueing Condition', 'Mean P(trans)', Titletext{mechanism});
            subplot(3,2,5);
            plot(mean(MMcwattraction,1));
            PostFigure([0.5, 3.5, 0, 1], 'Cueing Condition', 'Mean P(wheel)', Titletext{mechanism});
        end
        
        P.cuerate = CueRate; 

    end % mechanism
end  % task


for task = Tasks
    legendtext = {'With Mechanism', 'Without Any'};
    xlabel = {'Neu', 'Val', 'Inv'};
    PreFigure;
    index = 1;
    for mechanism = setdiff(Mechanisms, max(Mechanisms))
        if length(Mechanisms) < 4, subplot(1,2,index); end
        if length(Mechanisms) > 3 && length(Mechanisms) < 6, subplot(2,2,index); end
        if length(Mechanisms) > 5 && length(Mechanisms) < 8, subplot(2,3,index); end
        plotX = squeeze(meanAcc(task,mechanism,Cueconds));
        plot(Cueconds, plotX);
        hold on
        plotX = squeeze(meanAcc(task, max(Mechanisms), Cueconds));  % the last mechanism level is always: all mechanisms off
        plot(Cueconds, plotX, 'r');
        if (task==1), PostFigure([0.5, max(Cueconds)+0.5, 0, 90], 'Cueing Condition', 'Error (deg)', Titletext{mechanism}, legendtext); end
        if (task==2), PostFigure([0.5, max(Cueconds)+0.5, 0.5, 1], 'Cueing Condition', 'Accuracy', Titletext{mechanism}, legendtext); end
        xticklabels(xlabel(Cueconds));
        index = index + 1;
    end
end

out = 1; 


%
% %%% Save results
% if E.saveResults == 1
%     fid = fopen(['IMSim.SetsizeCueing', mat2str(PreRetro), '.dat'], 'w');
%     for id = 1:E.nsubj
%         for setsize = 1:E.maxsetsize
%             for cueing = 1:3
%                 fprintf(fid, '%d %d %d %d  ', id, setsize, cueing, Mdevobs(id, cueing, setsize));
%                 for ii = 1:length(indVar)
%                     fprintf(fid, '%d ', ParX(id, ii));
%                 end
%                 if (fitMM == 1), fprintf(fid, '%d %d %d %d', MMtranspos(id, cueing, setsize), MMguessing(id, cueing, setsize), MMcwattraction(id, cueing, setsize), MMSD(id, cueing, setsize)); end
%                 fprintf(fid, '\n');
%             end
%         end
%     end
%     fclose(fid);
% end




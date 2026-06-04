function [] = Refreshing(Model, setsize, fitMM, fitIMSim)
% Guided refreshing paradigm (Souza, Rerko, & Oberauer, 2015, ANYAS):
% Holding 6-item arrays in mind, with 4 successive cues guiding refreshing
% to individual items, so that 1 item is refreshed 2x, 2 items are refreshed 1x, and the remaining 3 items are refreshed 0x. 
% Memory is tested with continuous reproduction. 

global E
global C
global P

E.cuevalidity = 1/setsize;
option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

% Create refreshing sequences

abstractSeq = [1 2 1 3; 2 1 3 1];  % templates for the sequences of 3 items that receive refreshing cues
permutations = perms(1:setsize);  % all permutations of orders of the 6 items in the array
refreshed = [permutations(:, abstractSeq(1,:)); permutations(:, abstractSeq(2,:))]; % using abstract sequences as indices to pull out item numbers to be assigned to the template numbers
targetrefreshed = refreshed == 1; % item 1 is the target, so whenever item 1 is refreshed, the target is refreshed
Ntargetref = sum(targetrefreshed, 2);  % how many times is the target refreshed in each refreshing sequence? 
for refreshings = 0:2  % for each condition of "number of refreshings" (0, 1, and 2), ... 
    C.RefSequence(refreshings+1).seq = refreshed(Ntargetref == refreshings, :);  % ... pull out the subset of refreshing sequences in which the target is refreshed the right number of times
end

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff; 

% initialize parameters of mixture model, IM, observed mean Deviation

Mdevobs = zeros(E.nsubj, 3); % 0, 1, or 2 refreshings
MdevLastRef = NaN(E.nsubj, 2, 4); % for 1 or 2 refreshings: last refreshing in refreshing position
MMSD = zeros(E.nsubj, 3);     % SD parameter from Mixture Model
MMguessing = zeros(E.nsubj, 3);  % P(guessing) parameter from Mixture Model
MMtranspos = zeros(E.nsubj, 3);  % P(tranposition) parameter from Mixture Model

IMSimparms = zeros(E.nsubj, 6);   % IM parameters

Conditions = zeros(E.nsubj*3*E.ntrials, 1);   % matrix of conditions to be run
Array = zeros(E.nsubj*3*E.ntrials, setsize+1);  % arrays simulated
Target = zeros(1,E.nsubj*3*E.ntrials);          % Targets
Response = zeros(1,E.nsubj*3*E.ntrials);        % Responses
%Apred = zeros(E.nsubj*3*E.ntrials,360);        

[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x setsize x [1:360]

tcount = 1; %trial count

for id = 1:E.nsubj
    
        % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']); 
    end
    
    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);
    
    fdistance = zeros(1,E.ntrials);   % feature distance between target and response
    rt = zeros(1,E.ntrials);           % response times
    lastrefreshed = zeros(1,E.ntrials);   % last refreshed item
    Probedpos = zeros(E.ntrials,1);       % number of probed item (= probed position in the array)
    Pangle = zeros(E.ntrials,setsize);    % angle of the items' location in their circular arrangement in the array
    Cangle = zeros(E.ntrials,setsize+1);   % angles of the items' colors in the color wheel
    Targ = zeros(E.ntrials,1);             % target angles
    Resp = zeros(E.ntrials,1);             % response angles
    Setsize = zeros(E.ntrials,1);          % set size is important, too (except here it isn't because it is always 6)
    
    for refreshings = 1:4  % 1: no refreshing cues at all; 2: 0 refreshings, 3: 1 refreshing, 4: 2 refreshings
        
        for trial = 1:E.ntrials
            
            Conditions(tcount,1) = refreshings;
            cueing = 1 + 3*(refreshings>1);   %cueing = 4 -> refreshing
            
            output = Model(P, setsize, cueing, refreshings-1);  
            
            fdistance(trial) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
            rt(trial) = output.rt;
            lastrefreshed(trial) = output.lastrefreshed;
            Array(tcount,:) = [output.F(1:setsize), output.CWcolor];  %add the color in the color wheel closest to target location in column setsize+1                
            Target(tcount) = output.F(1);
            Response(tcount) = output.response;
            tcount = tcount+1;
            
            %collect data for further modeling
            Probedpos(trial) = output.L(1);
            Pangle(trial,:) = output.L(1:setsize);
            Cangle(trial,1:setsize) = output.F(1:setsize);
            Cangle(trial,setsize+1) = output.CWcolor;
            Targ(trial) = output.F(1);
            Resp(trial) = output.response;
            Setsize(trial) = setsize;
            
        end
        
        
        % aggregate data within each design cell (number of refreshings x position of last refreshing in the sequence) 
        Mdevobs(id, refreshings) = mean(abs(fdistance));  %mean deviation
        if refreshings>2
            for lastref = 1:4
                MdevLastRef(id, refreshings-2, lastref) = mean(abs(fdistance(lastrefreshed==lastref)));
            end
        end
        
        disp('    ID       RefCond  Trial      Error ');
        disp([id, refreshings, trial, Mdevobs(id, refreshings)]);
        
        condData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for mixture-model fitting: separate data structures for each condition
        condData.cueing = refreshings>2;
        condData.preretro = 2;
        condD.setsize = Setsize;
        condD.response = Resp;
        condD.L = round(C.Location(Pangle));
        condD.Color = Cangle;
        condD.cueing = refreshings>2;
        condD.preretro = 2;
        
        if refreshings == 1
            Data = condData;
            D = condD;
        else   % concatenate the data structures for individual conditions into a data structure that includes all conditions
            f = fieldnames(Data);
            for i = 1:length(f)
                Data.(f{i}) = [Data.(f{i}); condData.(f{i})];
            end
            ff = fieldnames(D);
            for i = 1:length(ff)
                D.(ff{i}) = [D.(ff{i}); condD.(ff{i})];
            end
        end
        
        % fit Mixture Model
        if fitMM
            startparms = [15, .1, .1];
            lb = [eps, 0, 0]; ub = [90, 1, 1];
            npar = 3;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture
            MMloglik = 500000;
            itercount = 0;
            while MMloglik > 400000
                [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, condData, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                itercount = itercount + 1;
            end
            disp('    ID        RefCond   Iter      MM-log(Lik)/1000 ');
            disp([id, refreshings, itercount, MMloglik./1000]);
            MMSD(id, refreshings) = MMparms(1);
            MMguessing(id, refreshings) = MMparms(2);
            if npar == 3, MMtranspos(id, refreshings) = MMparms(3); end
        end
        
    end % for cueing
    
    
    %fit IMSim
    if fitIMSim
        startparms = [0.05, 1.5, 3, 10, 20, 0.5];  %X, Y, s, kappa, kappafocus, Creduction
        npar = 6;
        lb = zeros(1,npar); ub = [5, 5, 20, 90, 90, 1];
        IMSimloglik = 500000;
        itercount = 0;
        while IMSimloglik > 400000
            [IMSimparms(id,:), IMSimloglik] = fminsearchbnd(@(x) IMFit(C.x, D, 2), startparms, lb, ub, option);
            itercount = itercount + 1;
        end
        disp([id, itercount, IMSimloglik]);
        pred1 = IMSim(IMSimparms(id,:), D, 1);
        Dev = abs(wrap(repmat(1:360, size(D.response,1), 1) - repmat(D.Color(:,1), 1, 360), 180));
        predDev1 = sum(Dev .* pred1, 2); % Weights deviation from each possible angle with that angle's predicted probability (of being the response). Computes weighted average deviation for predicted SD
        predMDevIMSim1(id, :) = aggregate(D.setsize, predDev1);
        
        pred2 = IMSim([X, Y, s, 2*kappa, 2*kappaf, r], D, 1);
        predDev2 = sum(Dev .* pred2, 2); % Weights deviation from each possible angle with that angle's predicted probability (of being the response). Computes weighted average deviation for predicted SD
        predMDevIMSim2(id, :) = aggregate(D.setsize, predDev2);
        
    end
    
end  % for ID

% Plot Mean(Deviation) as functions of number of refreshings, and last
% refreshing position

PreFigure
subplot(1,2,1);
plotvector = mean(Mdevobs(:, 2:4),1);
plot(0:2, plotvector);
hold on
plot([-0.2, 2.2], repmat(mean(Mdevobs(:,1)), 1, 2), '-b');  % no-refreshing-cue condition as blue line
PostFigure([-0.2, 2.2, 0, 80], 'N(Refreshings)', 'Deviation', 'Blue = No-Refreshing Cue');
subplot(1,2,2);
plotvector = squeeze(mean(MdevLastRef,1));
plot(1:4, plotvector');
hold on
plot(1:4, repmat(mean(Mdevobs(:,2)), 1, 4), '-r');  % zero refreshings baseline as red line
PostFigure([0.5, 4.5, 0, 80], 'Last Refreshing Position', 'Deviation', 'Red = 0 Refreshings', {'1 Ref', '2 Ref'});
hold off

% Plot response distributions
ResponseDistribC(Array, Target, Response, Conditions, {'No Ref', '0 Ref', '1 Ref', '2 Ref'});

% Plot Mixture Model Parameters over Conditions
if fitMM
    MMPm = 1 - MMtranspos - MMguessing;
    PreFigure
    subplot(2,2,1);
    plot(0:2, mean(MMSD(:, 2:4),1));
    hold on
    plot([-0.2, 2.2], repmat(mean(MMSD(:,1)), 1, 2), '-b');  % no-refreshing-cue condition as blue line
    PostFigure([-0.2, 2.2, 0, max(mean(MMSD,1))+0.5], 'N(Refreshings)', 'Mean SD', 'SD from Bays Mixture');
    subplot(2,2,2);
    plot(0:2, mean(MMPm(:, 2:4),1));
    hold on
    plot([-0.2, 2.2], repmat(mean(MMPm(:,1)), 1, 2), '-b');  % no-refreshing-cue condition as blue line
    PostFigure([-0.2, 2.2, 0, 1], 'N(Refreshings)', 'Mean P(m)', '"P(mem) Bays Mixture');
    subplot(2,2,3);
    plot(0:2, mean(MMguessing(:, 2:4),1));
    hold on
    plot([-0.2, 2.2], repmat(mean(MMguessing(:,1)), 1, 2), '-b');  % no-refreshing-cue condition as blue line
    PostFigure([-0.2, 2.2, 0, 0.5], 'N(Refreshings)', 'Mean P(guess)', '"P(guess) Bays Mixture');
    subplot(2,2,4);
    plot(0:2, mean(MMtranspos(:, 2:4),1));
    hold on
    plot([-0.2, 2.2], repmat(mean(MMtranspos(:,1)), 1, 2), '-b');  % no-refreshing-cue condition as blue line
    PostFigure([-0.2, 2.2, 0, 1], 'N(Refreshings)', 'Mean P(trans)', '"P(transpos) Bays Mixture');
end


if fitIMSim
    disp('      X       Y          s       kappa      kappaf      r');
    disp(mean(IMSimparms, 1));
    disp('      X       Y          s       kappa      kappaf      r');
    disp(std(IMSimparms, 1));
end


%%% Save results
fid = fopen('IMSim.Refreshing.dat', 'w');
for id = 1:E.nsubj
        for refreshings = 1:4
            fprintf(fid, '%d %d %d  ', id, refreshings, Mdevobs(id, refreshings));
            if (fitMM == 1), fprintf(fid, '%d %d %d', MMtranspos(id, refreshings), MMguessing(id, refreshings), MMSD(id, refreshings)); end
            fprintf(fid, '\n');
        end
end
fclose(fid); 

%%% Save results (with last-refreshed item)
fid = fopen('IMSim.RefreshingLast.dat', 'w');
for id = 1:E.nsubj
    for lastref = 1:4
        for refreshings = 1:2
            fprintf(fid, '%d %d %d %d  ', id, refreshings+2, lastref, MdevLastRef(id, refreshings, lastref));
            fprintf(fid, '\n');
        end
    end
end
fclose(fid); 

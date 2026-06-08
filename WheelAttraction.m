function D = WheelAttraction(Model, setsize, fitMM, fitIMSim)
% Simulates the color-wheel attraction effect (Experiment 6 in Souza,
% Rerko, & Oberauer, 2016, JEP:HPP): Arrays of 6 colors, continuous
% reproduction on a color wheel or a grey wheel (the color becomes visible
% only in the target object once the mouse is moved into the wheel).
% Retro-cue (= delay of onset of wheel) vs. no cue (= onset of wheel
% simultaneous with the cue identifying the target).

global C
global E
global P

option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);
E.cuevalidity = 1;
C.nloc = 6;  % only 6 equally spaced locations were used

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff; 

% initialize parameters of mixture model, IM, observed mean
% Deviation

Mdevobs = zeros(E.nsubj, 2, 2); % Mean deviation (target feature, response) for subjects, cueing condition (no/valid), and wheel type (color, grey)
MMSD = zeros(E.nsubj, 2, 2);    % SD parameter of mixture model
MMguessing = zeros(E.nsubj, 2, 2);  % P(guess) parameter of mixture model
MMtranspos = zeros(E.nsubj, 2, 2);  % P(transposition) parameter of mixture model
MMcwattraction = zeros(E.nsubj, 2, 2);  % P(color wheel attraction) parameter of mixture model
MRT = zeros(E.nsubj, 2, 2);   % subj, cueing condition, wheel condition

IMSimparms = zeros(E.nsubj, 6);

Conditions = zeros(E.nsubj*2*2*E.ntrials, 1);
Array = zeros(E.nsubj*2*2*E.ntrials, setsize+1);
Target = zeros(1,E.nsubj*2*2*E.ntrials);
Response = zeros(1,E.nsubj*2*2*E.ntrials);


[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x maxsetsize x [1:360]

tcount = 1; %trial count

for id = 1:E.nsubj
    
    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end
    
    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);
    
    for wheel = 1:2   % 1: color wheel, 2: grey wheel
        
        E.wheel = wheel; % depending on the wheel condition, set the color wheel encoded into WM at test to be present (1) or absent (2)
        fdistance = zeros(1,E.ntrials);
        rt = zeros(1,E.ntrials);
        Probedpos = zeros(E.ntrials,1);
        Pangle = zeros(E.ntrials,setsize);
        Cangle = zeros(E.ntrials,setsize);
        Targ = zeros(E.ntrials,1);
        Resp = zeros(E.ntrials,1);
        Setsize = zeros(E.ntrials,1);
        
        for cueing = 1:2  % 1: no cue (0 s delay), 2: valid cue (1 s delay),
            C.maxAdrift = 0; 
            for trial = 1:E.ntrials
                
                Conditions(tcount, 1) = (cueing-1)*2 + wheel;  % 1 = no-cue, colorwheel, 2 = no-cue, grey wheel, ...
                
                output = Model(P, setsize, cueing);
                
                Array(tcount,:) = [output.F(1,1:setsize), output.CWcolor];  %add the color in the color wheel closest to target location in column setsize+1
                fdistance(trial) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
                rt(trial) = output.rt;
                Target(tcount) = output.F(1,1);
                Response(tcount) = output.response;
                tcount = tcount+1;
                
                %collect data for further modeling
                Probedpos(trial) = output.L(1);
                Pangle(trial,:) = output.L(1:setsize);
                Cangle(trial,1:setsize) = output.F(1,1:setsize);
                Cangle(trial,setsize+1) = output.CWcolor;
                Targ(trial) = output.F(1,1);
                Resp(trial) = output.response;
                Setsize(trial) = setsize;
                
            end
            
            Mdevobs(id, cueing, wheel) = mean(abs(fdistance));  %mean deviation
            MRT(id, cueing, wheel) = mean(rt); 
            
            condData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for mixture-model fitting (separate structures for each condition)
            condData.cueing = repmat(cueing, length(Response), 1);
            condData.preretro = 2;
            condData.wheel = wheel;
            condD.setsize = Setsize;
            condD.response = Resp;
            condD.L = Pangle;
            condD.Color = Cangle;
            condD.cueing = repmat(cueing, length(Response), 1);
            condD.preretro = 2;
            condD.wheel = wheel;
            
            % Concatenate the condition-specific data to a data structure
            % across all conditions (for IM modelling)
            if wheel == 1 && cueing == 1
                D = condD;
            else   % concatenate the data structures
                ff = fieldnames(D);
                for i = 1:length(ff)
                    D.(ff{i}) = [D.(ff{i}); condD.(ff{i})];
                end
            end
            
            % fit Mixture Model
            if fitMM
                startparms = [15, .1, .1, .1];
                lb = [eps, 0, 0, 0]; ub = [90, 1, 1, 1];
                npar = 4;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture, 4 -> Souza & Oberauer mixture (including colorwheel attraction)
                MMloglik = 500000;
                itercount = 0;
                while MMloglik > 400000
                    [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, condData, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                    itercount = itercount + 1;
                end
                disp('    ID        Setsize    Cueing   iter      LL/1000   SD        Pmem      Ptrans    Pwheel ');
                disp([id, setsize, cueing, itercount, MMloglik/1000, MMparms(1), 1-sum(MMparms(2:4)), MMparms(3), MMparms(4)]);
                MMSD(id, cueing, wheel) = MMparms(1);
                MMguessing(id, cueing, wheel) = MMparms(2);
                if npar > 2, MMtranspos(id, cueing, wheel) = MMparms(3); end
                if npar > 3, MMcwattraction(id, cueing, wheel) = MMparms(4); end
            end
            
            disp('    ID        cueing    wheel    error ');
            disp([id, cueing, wheel, Mdevobs(id, cueing, wheel)]);

        end % for cueing
    end %for wheel
    
    
    
    %fit IMSim
    if fitIMSim
        startparms = [0.05, 1.5, 3, 10, 20, 0.5];  %X, Y, s, kappa, kappafocus, Creduction
        npar = 6;
        lb = zeros(1,npar); ub = [5, 5, 20, 90, 90, 1];
        IMSimloglik = 500000;
        itercount = 0;
        while IMSimloglik > 400000
            [IMSimparms(id,:), IMSimloglik] = fminsearchbnd(@(x) IMSim(x, D, 2), startparms, lb, ub, option);
            itercount = itercount + 1;
        end
        disp([id, itercount, IMSimloglik]);
        pred = IMSim(IMSimparms(id,:), D, 1);
        Dev = abs(wrap(repmat(1:360, size(D.response,1), 1) - repmat(D.Color(:,1), 1, 360), 180));
        predDev = sum(Dev .* pred, 2); % Weights deviation from each possible angle with that angle's predicted probability (of being the response). Computes weighted average deviation for predicted SD
        predMDevIMSim(id, :) = aggregate(D.setsize, predDev);
        
        pred2 = IMSim([X, Y, s, 2*kappa, 2*kappaf, r], D, 1);
        predDev2 = sum(Dev .* pred2, 2); % Weights deviation from each possible angle with that angle's predicted probability (of being the response). Computes weighted average deviation for predicted SD
        predMDevIMSim2(id, :) = aggregate(D.setsize, predDev2);
        
    end
    
end  % for ID

% Plot Mean(Deviation) as functions of Wheel
legendtext = {'Color Wheel', 'Grey Wheel'};

PreFigure
subplot(1,2,1);
plotvector = squeeze(mean(Mdevobs,1));
plot(plotvector);
PostFigure([0.8, 2.2, 0, 1.05*max(max(plotvector))], 'Cue Condition', 'Deviation', 'Mean Deviation', legendtext);
subplot(1,2,2);
plotvector = squeeze(mean(MRT,1));
plot(plotvector);
PostFigure([0.8, 2.2, 0, 1.05*max(max(plotvector))], 'Cue Condition', 'RT (ms)', 'Mean RT', legendtext);

D.Mdevobs = Mdevobs;
D.MRT = MRT;

% Plot response distributions
CondCodes = {'No Cue, Color Wheel', 'No Cue, Grey Wheel', 'Valid Cue, Color Wheel', 'Valid Cue, Grey Wheel'};
ResponseDistribC(Array, Target, Response, Conditions, CondCodes);

% Plot Mixture Model Parameters over Wheel condition
if fitMM
    MMPm = 1 - MMtranspos - MMguessing - MMcwattraction;
    Xticks = {'Color', ' ', 'Grey'};  % Matlab insists on placing a tick every 0.5 steps, so have to insert a blank tick
    PreFigure
    subplot(3,2,1);
    plot(squeeze(mean(MMSD,1))');
    set(gca,'XtickLabel', Xticks);
    PostFigure([0.8,2.2, 0, max(max(mean(MMSD,1)))+0.5], 'Wheel Condition', 'Mean SD', 'SD from Bays Mixture', {'No Cue', 'Retro-Cue'});
    subplot(3,2,2);
    plot(squeeze(mean(MMPm,1))');
    set(gca,'XtickLabel', Xticks);
    PostFigure([0.8,2.2, 0, 1], 'Wheel Condition', 'Mean P(m)', '"P(mem) Bays Mixture');
    subplot(3,2,3);
    plot(squeeze(mean(MMguessing,1))');
    set(gca,'XtickLabel', Xticks);
    PostFigure([0.8,2.2, 0, 0.5], 'Wheel Condition', 'Mean P(guess)', 'P(guess) Bays Mixture');
    subplot(3,2,4);
    plot(squeeze(mean(MMtranspos,1))');
    set(gca,'XtickLabel', Xticks);
    PostFigure([0.8,2.2, 0, 0.5], 'Wheel Condition', 'Mean P(trans)', 'P(transpos) Bays Mixture');
    subplot(3,2,5);
    plot(squeeze(mean(MMcwattraction,1))');
    PostFigure([0.8,2.2, 0, 0.5], 'Wheel Condition', 'Mean P(CW)', 'P(CW) Bays Mixture');
    
    % Figure layout as in Souza et al (2016)
    meanPm = squeeze(mean(MMPm, 1));
    meanCWA = squeeze(mean(MMcwattraction, 1));
    Xticks = {'No Cue', ' ', 'Cue'};
    PreFigure;
    subplot(2,2,1);
    plot(1:2, meanPm(:,1)); % color-wheel
    set(gca,'XtickLabel', Xticks);
    PostFigure([0.8, 2.2, 0, 1], [], 'Mean P(m)', 'Color Wheel');
    subplot(2,2,2);
    plot(1:2, meanPm(:,2));
    set(gca,'XtickLabel', Xticks);
    PostFigure([0.8,2.2, 0, 1], [], 'Mean P(m)', 'Grey Wheel');
    subplot(2,2,3);
    plot(1:2, meanCWA(:,1)); % color-wheel
    set(gca,'XtickLabel', Xticks);
    PostFigure([0.8,2.2, -0.03, 0.5], [], 'Mean P(CW)', 'Color Wheel');
    subplot(2,2,4);
    plot(1:2, meanCWA(:,2));
    set(gca,'XtickLabel', Xticks);
    PostFigure([0.8,2.2, -0.03, 0.5], [], 'Mean P(CW)', 'Grey Wheel');
    
    D.MMSD = MMSD;
    D.MMpm = MMPm;
    D.MMguessing = MMguessing;
    D.MMtranspos = MMtranspos;
    D.MMcwattraction = MMcwattraction;

end

if fitIMSim
    disp('      X       Y          s       kappa      kappaf      r');
    disp(mean(IMSimparms, 1));
    disp('      X       Y          s       kappa      kappaf      r');
    disp(std(IMSimparms, 1));
end



%%% Save results
if E.saveResults == 1
    fid = fopen('IMSim.WheelAttraction.dat', 'w');
    for id = 1:E.nsubj
        for wheel = 1:2
            for cueing = 1:2
                fprintf(fid, '%d %d %d %d  ', id, wheel, cueing, Mdevobs(id, cueing, wheel));
                if (fitMM == 1), fprintf(fid, '%d %d %d %d', MMtranspos(id, cueing, wheel), MMguessing(id, cueing, wheel), MMcwattraction(id, cueing, wheel), MMSD(id, cueing, wheel)); end
                fprintf(fid, '\n');
            end
        end
    end
    fclose(fid);
    fid2 = fopen('IMSim.WheelAttractionDetails.dat', 'w');
    for trial = 1:size(Conditions,1)
        fprintf(fid2, '%d  %d %d ', Conditions(trial), Target(trial), Response(trial));
        for item = 1:(setsize+1)   % column setsize+1 contains color spatially closest to the target in the color wheel
            fprintf(fid2, '%d ', Array(trial,item));
        end
        fprintf(fid2, '\n');
    end
    fclose(fid2);
end



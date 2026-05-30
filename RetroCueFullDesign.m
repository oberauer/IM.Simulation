function [] = RetroCueFullDesign(Model, indVar, maxIndVar, fitMM)
% Simulation of retro-cue effect varying task (continuous reproduction vs. change detection, cue validity, and the
% presence/absence of the 3 mechanisms (strengthening, removal, and perceptual interference

global P
global E
global C

E.PreRetro = 2;
setsize = 6;
E.cuevalidity = 1;
E.material = 2;

IMprepareRecog; % set up criterion for expected size of change 
Ptype = [1 1 2 3];  % 2 x positive, 1 x new, 1 x intrusion
option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

cueingStrength = P.cueingStrength;
testdisplay = P.testdisplay;
removalTau = P.removalTau;
removalGain = P.removalGain;

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff; 

for task = 1:2

    E.test = task;
    E.material = task; % for change detection, use small set of highly distinct stimuli

    for cueTarget = 1:2  % 1 = last cued item, 2 = next-to-last cued item

        meanAcc = zeros(2,2,2,2);

        for str = 1:2

            for rem = 1:2

                for inter = 1:2

                    Mdevobs = NaN(E.nsubj, 2);  % mean observed deviation of responses from target feature for subjects, cueing conditions,
                    Mrt = NaN(E.nsubj, 2);      % mean RT
                    MMSD = NaN(E.nsubj, 2);     % Mixture Model SD parameter
                    MMguessing = NaN(E.nsubj, 2);  % Mixture Model Guessing parameter
                    MMtranspos = NaN(E.nsubj, 2);  % Mixture Model Transposition parameter (swap error proportion)
                    MMcwattraction = NaN(E.nsubj, 2);  % Mixture Model colorwheel-attraction strength parameter
                    Pcorrect = NaN(E.nsubj, 2, 3);       % mean accuracy in change detection for 3 cueing conditions and 3 probe types

                    % pre-define container matrices for memory arrays, target features, and
                    % responses
                    Array = zeros(E.nsubj*3*E.ntrials, setsize);
                    Target = zeros(1,E.nsubj*3*E.ntrials);
                    Response = zeros(1,E.nsubj*3*E.ntrials);

                    % for each set-size level, for each trial, generate a vector of 360 color
                    % values coding the colors on the wheel
                    [aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x E.maxsetsize x [1:360]

                    tcount = 1; %trial count

                    for id = 1:E.nsubj

                        % extract parameter values for each subject - for those parameters that vary between subjects
                        for ii = 1:length(indVar)
                            eval(['P.', indVar{ii}, ' = ParX(id, ii);']);
                        end
                        if str == 1, P.cueingStrength = cueingStrength; else, P.cueingStrength = 0; end
                        if rem == 1, P.removalTau = 0; P.removalGain = -100; else, P.removalTau = removalTau; P.removalGain = removalGain; end
                        if inter == 1, P.testdisplay = 0; else, P.testdisplay = testdisplay; end

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

                        for cueingIdx = 1:2  % neutral, valid (single cue or AB cue)

                            if cueingIdx == 1, cueing = 1; end
                            if cueingIdx == 2, cueing = cueingIdx + 2*(cueTarget==2); end

                            for trial = 1:E.ntrials

                                if task == 2
                                    E.ptype = Ptype(mod(trial,4)+1);
                                end
                                if cueTarget == 2, C.RefSequence(1).seq = [1, 1+randperm(setsize-1, 1)]; end

                                output = Model(P, setsize, cueing);   % here the model is run!

                                Array(tcount,:) = output.F(1,1:setsize);  % record of the array on this trial
                                Target(tcount) = output.F(1,1);    % target feature

                                if task == 1
                                    Response(tcount) = output.response;   % response feature
                                    fdistance(trial) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
                                end
                                if task == 2
                                    Response(tcount) = output.response(1,:);  % the first entry of response is the actual response
                                    Probetype(tcount) = E.ptype;
                                    %delta = output.response(2);
                                end
                                rt(trial) = output.rt;            % response time

                                tcount = tcount+1;   % trial counter is incremented

                                if task == 1
                                    %collect data for further modeling with Mixture Model
                                    Probedpos(trial) = output.L(1);    % probed position
                                    Pangle(trial,:) = output.L(1:setsize);  % spatial position angles
                                    Cangle(trial,1:setsize) = output.F(1,1:setsize);  % array feature angles
                                    Cangle(trial,setsize+1) = output.CWcolor;  % color-wheel feature angles
                                    Targ(trial) = output.F(1,1);           % target feature
                                    Resp(trial) = output.response;       % response feature
                                    Setsize(trial) = setsize;
                                end

                            end

                            if task == 1
                                Mdevobs(id, cueingIdx) = mean(abs(fdistance));  %mean deviation (averaged over trials)
                            end
                            if task == 2
                                Pyes = 2-Response;  % Yes/No: response = 1/2
                                Pcorrect(id, cueingIdx, 1) = mean(Pyes(Probetype==1));
                                for ptype = 2:3
                                    Pcorrect(id, cueingIdx, ptype) = mean(1-Pyes(Probetype==ptype));
                                end
                            end
                            Mrt(id, cueingIdx) = mean(rt);

                            if fitMM == 1 && task == 1
                                ssData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting
                                %             ssData.cueing = repmat(cueing, length(Response), 1);
                                %             ssData.preretro = 2;

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
                                disp('    id        cueing    strengthening  removal   vis. interference  Deviance     SD(deg)   P(guess)  P(trans)   P(CW)');
                                disp([id, cueingIdx, str, rem, inter, MMloglik/1000, MMparms]);
                                MMSD(id, cueingIdx) = MMparms(1);
                                MMguessing(id, cueingIdx) = MMparms(2);
                                if npar > 2, MMtranspos(id, cueingIdx) = MMparms(3); end
                                if npar > 3, MMcwattraction(id, cueingIdx) = MMparms(4); end
                            else
                                disp('    id    task #cues cueing strength rem inter ');
                                disp([id, task, cueTarget, cueingIdx, str, rem, inter]);
                            end

                        end % for cueing

                    end  % for ID


                    if task == 1, Accuracy = Mdevobs; end
                    if task == 2, Accuracy = squeeze(mean(Pcorrect,3)); end  % average over probe types
                    
                    % PreFigure;
                    % plotvector = mean(Accuracy,1);
                    % plot(plotvector);
                    % titletext = ['Cueval = ', mat2str(round(E.cuevalidity,2)), '; Str = ', mat2str(consolid), '; Rem = ', mat2str(rem), '; Inter = ', mat2str(inter)];
                    % PostFigure([0.5, 3.5, 0, 1.05*max(max(Accuracy))], 'Cueing Condition', 'Accuracy', titletext);

                    meanAcc(str, inter, rem, :) = mean(Accuracy, 1);

                    % Plot Mixture Model Parameters over Setsize
                    if fitMM && task == 1
                        MMPm = 1 - MMtranspos - MMguessing - MMcwattraction;
                        PreFigure
                        subplot(3,2,1);
                        plot(mean(MMSD,1));
                        PostFigure([0.5, 3.5, max(max(mean(MMSD,1)))+0.5], 'Cueing Condition', 'Mean SD', 'SD from Bays Mixture');
                        subplot(3,2,2);
                        plot(mean(MMPm,1));
                        PostFigure([0.5, 3.5, 0, 1], 'Cueing Condition', 'Mean P(m)', 'P(mem)');
                        subplot(3,2,3);
                        plot(mean(MMguessing,1));
                        PostFigure([0.5, 3.5, 0, 0.5], 'Cueing Condition', 'Mean P(guess)', 'P(guess)');
                        subplot(3,2,4);
                        plot(mean(MMtranspos,1));
                        PostFigure([0.5, 3.5, 0, 1], 'Cueing Condition', 'Mean P(trans)', 'P(transpos)');
                        subplot(3,2,5);
                        plot(mean(MMcwattraction,1));
                        PostFigure([0.5, 3.5, 0, 1], 'Cueing Condition', 'Mean P(wheel)', 'P(wheel attraction)');
                    end
                end % interference
            end  % removal
        end  % consolidation

        legendtext = {'No Str., no Inter', 'Str., no Inter' ,'No Str., Inter', 'Str., Inter'};
        CueingSchema = {'Target', 'Target-Other'}; 
        PreFigure;

        subplot(1,2,1);
        plotX = squeeze(meanAcc(1,1,1,:));  % consolid off, interference off, removal off
        plot(1:2, plotX);
        hold on
        plotX = squeeze(meanAcc(2,1,1,:));  % consolid on, interference off, removal off
        plot(1:2, plotX, 'r');
        plotX = squeeze(meanAcc(1,2,1,:));  % consolid off, interference on, removal off
        plot(1:2, plotX, 'b');
        plotX = squeeze(meanAcc(2,2,1,:));  % consolid on, interference on, removal off
        plot(1:2, plotX, 'g');
        PostFigure([0.5, 2.5, 0, 1.1*max(max(plotX))], 'Cueing Condition', 'Accuracy', ['Cueing Schema = ', CueingSchema{cueTarget}, '; No Removal'], legendtext);
        
        subplot(1,2,2);
        plotX = squeeze(meanAcc(1,1,2,:));  % consolid off, interference off, removal on
        plot(1:2, plotX);
        hold on
        plotX = squeeze(meanAcc(2,1,2,:));  % consolid on, interference off, removal on
        plot(1:2, plotX, 'r');
        plotX = squeeze(meanAcc(1,2,2,:));  % consolid off, interference on, removal on
        plot(1:2, plotX, 'b');
        plotX = squeeze(meanAcc(2,2,2,:));  % consolid on, interference on, removal on
        plot(1:2, plotX, 'g');
        PostFigure([0.5, 2.5, 0, 1.1*max(max(plotX))], 'Cueing Condition', 'Accuracy', ['Cue Schema = ', CueingSchema{cueTarget}, '; Removal'], legendtext);

    end  % cueing schema (A vs. AB)

end  % task

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




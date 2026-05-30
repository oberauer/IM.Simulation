function [] = Consolidation(Model, fitMM)
% Time course of consolidation and masking effects: Ricker & Sandry (2018)

global P
global E
global C

if ~exist('fitMM'), fitMM = 0; end

setsize = 3;
C.tstep = 0.02;  % needs a finer time step for the fine temporal resolution of consolidation times and masks
E.wheel = 2;  % no color wheel
E.presentation = 2;   % sequential presentation
E.outsize = setsize;  % test all items!
E.forwardrecall = 1;  % forward order of recall
E.RI = 0;  % no additional retention interval after the last presentation interval
%if ~ismember('cRate', C.indVar), P.cRate = P.cRate * P.cRateFactor; end
P.cRate = P.cRate * P.cRateFactor; 
%P.encrate = P.encrate * P.encrateFactor; 
%P.closeFXrate = P.closeFXrate * P.closeFXrateFactor; 
Ptime = [0.2, 0.3, 0.4, 0.5, 0.7, 1.0];  % presentation + consolidation time (= SOA between items), including the 0.067 s of actual on-time in Ricker & Sandry (2018)
MaskSOA = [0.08, 0.17, 10];    % 10 for the no-mask condition in Ricker & Sandry.

option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

% initialize parameters of mixture model, IM, observed mean Deviation

Mdevobs = NaN(E.nsubj, length(MaskSOA), length(Ptime), setsize);  % id, mask SOA, presentation time, serial position
mStrength = NaN(E.nsubj, length(MaskSOA), length(Ptime), setsize);  % id, mask SOA, presentation time, serial position
MMSD = NaN(E.nsubj, length(MaskSOA), length(Ptime), setsize);  % id, mask SOA, presentation time, serial position
MMguessing = NaN(E.nsubj, length(MaskSOA), length(Ptime), setsize);  % id, mask SOA, presentation time, serial position
MMtranspos = NaN(E.nsubj, length(MaskSOA), length(Ptime), setsize);  % id, mask SOA, presentation time, serial position
MMcw = NaN(E.nsubj, length(MaskSOA), length(Ptime), setsize);  % id, mask SOA, presentation time, serial position

%For plotting response distributions:
Array = zeros(E.nsubj*length(Ptime)*length(MaskSOA)*E.ntrials, setsize);
Target = zeros(1,E.nsubj*length(Ptime)*length(MaskSOA)*E.ntrials);
Response = zeros(1,E.nsubj*length(Ptime)*length(MaskSOA)*E.ntrials);
Conditions = zeros(E.nsubj*length(Ptime)*length(MaskSOA)*E.ntrials, 1);

[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x setsize x [1:360]

tcount = 1; %trial count

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    if ismember('cRate', C.indVar), P.cRate = P.cRate * P.cRateFactor; end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;   % 8 dots on a ring created the mask
    CreateMapping(E.calibrateAmp==2);

    for mtime = 1:length(MaskSOA)
        for ptime = 1:length(Ptime)

            E.prestime = Ptime(ptime);
            E.MaskSOA = MaskSOA(mtime);
            fdistance = zeros(E.ntrials, setsize);  % distance (target, response)
            strength = zeros(E.ntrials, setsize);

            Probedpos = zeros(E.ntrials,1);
            Pangle = zeros(E.ntrials,setsize);  % spatial angles of item positions in the array
            Cangle = zeros(E.ntrials,setsize+1); % color angles in the color wheel
            Targ = zeros(E.ntrials,1);
            Resp = zeros(E.ntrials,1);
            Setsize = zeros(E.ntrials,1);

            for trial = 1:E.ntrials

                output = Model(P, setsize, 1);  % cueing = 1 (no cue)
                %fdistance(trial, :) = wrap(output.response-output.F(1,1:setsize), 180);   %calculate distance between response and true feature in feature space (degrees!)
                
                for outpos = 1:setsize
                    fdistance(trial, output.Inpos(outpos)) = wrap(output.response(outpos)-output.F(1,outpos), 180);   %calculate distance between response and true feature in feature space (degrees!)
                    RT(trial, output.Inpos(outpos)) = output.rt(outpos);
                end
                               
                strength(trial, :) = output.Strength;

                %collect data for further modeling - only the first item tested
                Probedpos(trial) = output.L(1);    % probed position
                Pangle(trial,:) = output.L(1:setsize);  % spatial position angles
                Cangle(trial,1:setsize) = output.F(1,1:setsize);  % array feature angles
                Cangle(trial,setsize+1) = output.CWcolor;  % color-wheel feature angles
                Targ(trial) = output.F(1, 1);
                Resp(trial) = output.response(1);
                Setsize(trial) = setsize;

                % data for plotting of error distribution
                Conditions(tcount, 1) = (mtime-1)*length(Ptime) + ptime;
                Array(tcount,:) = output.F(1,1:setsize);  % record of the array on this trial
                Target(tcount) = output.F(1,1);    % target feature
                Response(tcount) = output.response(1);   % response feature
                tcount = tcount+1;

            end

            Mdevobs(id, mtime, ptime, :) = mean(abs(fdistance), 1);  %mean deviation, average over trials
            mStrength(id, mtime, ptime, :) = mean(strength, 1);

            ssData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting
            % ssD.setsize = Setsize;
            % ssD.response = Resp;
            % ssD.Color = Cangle;

            % fit Mixture Model
            if fitMM
                startparms = [15, .1, .1, .1];
                lb = [eps, 0, 0, 0]; ub = [90, 1, 1, 1];
                npar = 4;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture, 4 -> Souza & Oberauer mixture (iC.ncluding color-wheel attraction)
                MMloglik = 500000;
                itercount = 0;
                while MMloglik > 400000
                    [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, ssData, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                    itercount = itercount + 1;
                end
                disp('    ID        MaskSOA   ConsTime  Iter      LL/1000   SD        Pmem');
                disp([id, mtime, ptime, itercount, MMloglik/1000, MMparms(1), 1-sum(MMparms(2:end))]);
                MMSD(id, mtime, ptime) = MMparms(1);
                MMguessing(id, mtime, ptime) = MMparms(2);
                if npar > 2, MMtranspos(id, mtime, ptime) = MMparms(3); end
                if npar > 3, MMcw(id, mtime, ptime) = MMparms(4); end
            else
                disp('    ID        MaskSOA   ConsTime   Error');
                disp([id, mtime, ptime, mean(Mdevobs(id, mtime, ptime, :))]);
            end


        end  % for ptime
    end %for mtime

end  % for ID

% Plot Mean(Deviation) as function of serial position, consolidation time, and mask SOA

% Errors
legendtext = vec2legend(1:setsize);
PreFigure;
for mtime = 1:length(MaskSOA)
    subplot(2,2,mtime);
    plotvector = squeeze(mean(Mdevobs(:,mtime,:,:),1))';
    plot(Ptime, plotvector);
    PostFigure([0, 1.1*max(Ptime), 0, 90], 'Consolidation Time', 'Deviation (Deg)', ['SOA = ', mat2str(MaskSOA(mtime))], legendtext);
end

% Errors as a function of consolidation time (panels) and mask-onset time
% (x-axis)
legendtext = vec2legend(1:setsize);
xTicks = vec2legend(MaskSOA);
xTicks{end} = 'None'; 
PreFigure;
for ptime = 1:length(Ptime)
    if length(Ptime) < 5, subplot(2,2,ptime); end
    if length(Ptime) > 4, subplot(3,2,ptime); end
    plotvector = squeeze(mean(Mdevobs(:,:,ptime,:),1));
    plot(1:length(MaskSOA), plotvector);
    PostFigure([0, 1.2*length(MaskSOA), 0, 90], 'Array-Mask SOA', 'Deviation (Deg)', ['Cons. Time = ', mat2str(Ptime(ptime))], legendtext);
    xticks(1:length(MaskSOA)); 
    xticklabels(xTicks);
end



%%% Save results
if E.saveResults == 1
    fid = fopen(['IMSim.Consolidation.dat'], 'w');
    for id = 1:E.nsubj
        for mtime = 1:length(MaskSOA)
            for ptime = 1:length(Ptime)
                for outpos = 1:setsize
                    fprintf(fid, '%d %d %d %d  %d  ', id, mtime, ptime, outpos, Mdevobs(id, mtime, ptime, outpos));
                    fprintf(fid, '\n');
                end
            end
        end
    end
    fclose(fid);
end



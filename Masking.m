function [] = Masking(Model, fitMM)
% Simulation of structure-mask experiment of Agaoglu et al. (2015)

global P
global E
global C

fitIM = 0;

E.presentation = 1;  % simultaneous
E.prestime = 0.01;   % 10 ms
E.RI = 1.5; 
SOA = [-0.1, -0.05, -0.01, 0, 0.02, 0.04, 0.05, 0.06, 0.08, 0.11, 0.15, 0.2]; % from Exp. 3 of Agaoglu et al. 

option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

% initialize parameters of mixture model, IM, observed mean Deviation
Mdevobs = NaN(E.nsubj, length(SOA));  % id, soa
CircSD = NaN(E.nsubj, length(SOA));  % id, soa
MMSD = NaN(E.nsubj, length(SOA));  % 
MMguessing = NaN(E.nsubj, length(SOA));  %
MMtranspos = NaN(E.nsubj, length(SOA));  % id, setsize, inpos, outpos

IMparms = zeros(E.nsubj, 6);

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli(3);
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, 4), 1:360);  %Colors = E.ntrials x set-size x [1:360]

tcount = 1; %trial count

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    maskFeatures = CreateStimuli(3);
    CreateMapping(E.calibrateAmp==2);

    for soaIdx = 1:length(SOA)
        
        E.MaskSOA = SOA(soaIdx); 
        fdistance = NaN(E.ntrials, 1);  % distance (target, response) 
        Probedpos = zeros(E.ntrials,1);
        Pangle = zeros(E.ntrials, 4);  % spatial angles of item positions in the array - the same for target and the three mask features
        Cangle = zeros(E.ntrials, 5); % feature angles of the target (1st position) and the 3 mask features (2:4) + one for the colorwheel, which is not used here
        Targ = zeros(E.ntrials,1);
        Resp = zeros(E.ntrials,1);
        Setsize = zeros(E.ntrials,1);

        for trial = 1:E.ntrials

            output = Model(P, 1, 1);  % setsize = 1,  cueing = 1 (no cue)
            fdistance(trial, 1) = wrap(output.response(1)-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)

            %collect data for further modeling - only the first item tested
            Probedpos(trial) = output.L(1);
            Pangle(trial,:) = output.L(1:1);
            Cangle(trial,1) = output.F(1,1);
            Cangle(trial,2:4) = maskFeatures;
            Targ(trial) = output.F(1);
            Resp(trial) = output.response(1);
            Setsize(trial) = 4;  %including the mask features

        end

         Mdevobs(id, soaIdx) = mean(abs(fdistance));  %mean deviation
         CircSD(id, soaIdx) = circ_std(deg2rad(fdistance));
            
            % fit Mixture Model
            if fitMM
                ssData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting
                startparms = [15, .1, .1, .1];
                lb = [eps, 0, 0, 0]; ub = [90, 1, 1, 1];
                npar = 3;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture, 4 -> Souza & Oberauer mixture (iC.ncluding color-wheel attraction)
                MMloglik = 500000;
                itercount = 0;
                while MMloglik > 400000
                    [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, ssData, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                    itercount = itercount + 1;
                end
                disp('    ID        MaskSOA   Iter      LL/1000   SD        P(unif.)     P(swap)');
                disp([id, SOA(soaIdx), itercount, MMloglik/1000, MMparms]);
                MMSD(id, soaIdx) = MMparms(1);
                MMguessing(id, soaIdx) = MMparms(2);
                MMtranspos(id, soaIdx) = MMparms(3); 
            else
                disp('    ID        MaskSOA   Error');
                disp([id, SOA(soaIdx), mean(Mdevobs(id, soaIdx))]);
            end

    end %for SOA

    %fit IM
    if fitIM
        startparms = [0.5, 0.5, 2, 10, 20, 0.5];  %B, A, s, P.kappa, P.kappafocus, Creduction
        npar = 6;
        lb = zeros(1,npar); ub = [5, 5, 20, 90, 90, 1];
        IMloglik = 500000;
        itercount = 0;
        while IMloglik > 400000
            [IMparms(id,:), IMloglik] = fminsearchbnd(@(x) IM(x, Data, 2), startparms, lb, ub, option);
            itercount = itercount + 1;
        end
        disp([id, itercount, IMloglik]);
        pred = IM(IMparms(id,:), Data, 1);
        Dev = abs(wrap(repmat(1:360, size(Data.response,1), 1) - repmat(Data.response, 1, 360), 180));
        predDev = sum(Dev .* pred, 2); % Weights deviation from each possible angle with that angle's predicted probability (of being the response). Computes weighted average deviation for predicted SD
        predMDevIM(id, :) = aggregate(Data.setsize, predDev);
    end


end  % for ID

% Plot performance as function of SOA

TransformedPerformance = 1 - Mdevobs./90; 

PreFigure;
subplot(2,2,1);
plotvector = mean(TransformedPerformance);  % average over subjects
plot(SOA, plotvector);
PostFigure([min(SOA), max(SOA), 0, 1], 'SOA', 'Transf. Perf.');

subplot(2,2,3);
plotvector = mean(MMguessing);  % average over subjects
plot(SOA, plotvector);
PostFigure([min(SOA), max(SOA), 0, 1], 'SOA', 'P(Uniform)');

subplot(2,2,4);
plotvector = mean(MMtranspos);  % average over subjects
plot(SOA, plotvector);
PostFigure([min(SOA), max(SOA), 0, 1], 'SOA', 'P(Swap)');

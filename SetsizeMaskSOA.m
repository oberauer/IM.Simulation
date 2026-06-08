function D = SetsizeMaskSOA(Model, fitMM)
% Simulation of Set-size and array-mask SOA (Bays et al., JoV 2011)
global P
global E
global C

option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);
fitIM = 0;
E.maxsetsize = 6; 
E.targetDim = 2;  % target dimension is orientation
E.material = 3; 
E.mask = 2;    % Bays et al. mask covering the entire field
SOA = [0.025, 0.05, 0.075, 0.1, 0.125, 0.3, 0.5, 1, 1.5, 2];

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff; 

% initialize parameters of mixture model, IM, observed mean Deviation

Mdevobs = NaN(E.nsubj, length(SOA), E.maxsetsize);  % mean observed deviation of responses from target feature for subjects, cueing conditions, and set size
Mprecision = NaN(E.nsubj, length(SOA), E.maxsetsize);  % mean precision (= 1/circ_std) of response deviations from target
Mrt = NaN(E.nsubj, length(SOA), E.maxsetsize);      % mean RT
MMSD = NaN(E.nsubj, length(SOA), E.maxsetsize);     % Mixture Model SD parameter
MMpm = NaN(E.nsubj, length(SOA), E.maxsetsize);     % Mixture Model Pm parameter
MMguessing = NaN(E.nsubj, length(SOA), E.maxsetsize);  % Mixture Model Guessing parameter
MMtranspos = NaN(E.nsubj, length(SOA), E.maxsetsize);  % Mixture Model Transposition parameter (swap error proportion)
MMcwattraction = NaN(E.nsubj, length(SOA), E.maxsetsize);  % Mixture Model colorwheel-attraction strength parameter
Mwact = NaN(E.nsubj, length(SOA), E.maxsetsize);       % mean activation in the binding weight matrix
ChancePrecision = 1./circ_std((pi*(-90:90)./180)');

IMparms = zeros(E.nsubj, 6);     % Parameters for Interference Mdel

% pre-define container matrices for memory arrays, target features, and
% responses
Array = zeros(E.nsubj*3*E.maxsetsize*E.ntrials, E.maxsetsize);
Target = zeros(1,E.nsubj*3*E.maxsetsize*E.ntrials);
Response = zeros(1,E.nsubj*3*E.maxsetsize*E.ntrials);

% for each set-size level, for each trial, generate a vector of 360 color
% values coding the colors on the wheel
[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, E.maxsetsize), 1:360);  %Colors = E.ntrials x E.maxsetsize x [1:360]

tcount = 1; %trial count

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']); 
    end
    
    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;   % a large number of orientations were used to create the mask in Bays et al.
    CreateMapping(E.calibrateAmp==2);
    
    
    for setsize = 1:E.maxsetsize
        
        % Initialize container vectors
        fdistance = zeros(1,E.ntrials);  % feature distance between response and target
        rt = zeros(1,E.ntrials);         % response time
        wact = zeros(1,E.ntrials);       % activation (summed strength of bindings) of binding weight matrix
        Probedpos = zeros(E.ntrials,1);  % Number of the tested (probed) spatial position
        Pangle = zeros(E.ntrials,E.maxsetsize);  % spatial angles of item positions in the array
        Cangle = zeros(E.ntrials,E.maxsetsize+1); % color angles in the color wheel
        Targ = zeros(E.ntrials,1);                % Target
        Resp = zeros(E.ntrials,1);                % Response
        Setsize = zeros(E.ntrials,1);             % Set size
        
        for soa = 1:length(SOA)  %
            
            E.MaskSOA = SOA(soa);
            E.RI = E.MaskSOA + 1;     
            
            for trial = 1:E.ntrials
                
                output = Model(P, setsize, 1);   % here the model is run!
                
                Array(tcount,:) = [output.F(1,1:setsize), zeros(1, E.maxsetsize-setsize)];  % record of the array on this trial (fill up the rest until max(setsize) with zeros
                Target(tcount) = output.F(1,1);    % target feature
                Response(tcount) = output.response;   % response feature
                rt(trial) = output.rt;            % response time
                if E.material == 3
                    fdist(1) = wrap(output.response-output.F(1), 180);   %calculate distance between response and first orientation of a bar
                    fdist(2) = wrap(output.response-(output.F(1)+180), 180);   %calculate distance between response and second (opposite) orientation of the bar
                    fdistance(trial) = fdist(find(abs(fdist) == min(abs(fdist)),1)); 
                else
                    fdistance(trial) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
                end
                wact(trial) = sum(sum(output.g)); % sum of activation in gating layer -> CDA?
                tcount = tcount+1;   % trial counter is incremented
                
                %collect data for further modeling with Mixture Model or IM
                Probedpos(trial) = output.L(1);    % probed position
                Pangle(trial,:) = output.L(1:E.maxsetsize);  % spatial position angles
                Cangle(trial,1:setsize) = output.F(1,1:setsize);  % array feature angles
                Cangle(trial,E.maxsetsize+1) = output.CWcolor;  % color-wheel feature angles
                Targ(trial) = output.F(1,1);           % target feature
                Resp(trial) = output.response;       % response feature
                Setsize(trial) = setsize;            % set size
                
            end
            
            Mdevobs(id, soa, setsize) = mean(abs(fdistance));  %mean deviation (averaged over trials)
            Mprecision(id, soa, setsize) = 1./circ_std((pi*fdistance./180)') - ChancePrecision; 
            Mrt(id, soa, setsize) = mean(rt);
            Mwact(id, soa, setsize) = mean(wact);
            
            ssData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting
            ssData.cueing = repmat(soa, length(Resp), 1);
            ssData.preretro = E.PreRetro;
            ssD.setsize = Setsize;
            ssD.response = Resp;
            ssD.L = round(C.Location(Pangle));
            ssD.Color = Cangle;
            ssD.cueing = repmat(soa, length(Resp), 1);
            ssD.preretro = E.PreRetro;
            
            if setsize == 1
                Data = ssData;
                D = ssD;
            else   % concatenate the data structures
                f = fieldnames(Data);
                for i = 1:length(f)
                    Data.(f{i}) = [Data.(f{i}); ssData.(f{i})];
                end
                ff = fieldnames(D);
                for i = 1:length(ff)
                    D.(ff{i}) = [D.(ff{i}); ssD.(ff{i})];
                end
            end
            
            % fit Mixture Model
            if fitMM
                startparms = [15, .1, .1, .1];
                lb = [eps, 0, 0, 0]; ub = [90, 1, 1, 1];
                npar = 3;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture, 4 -> Souza & Oberauer mixture (iC.ncluding color-wheel attraction)
                MMloglik = 500000;
                itercount = 0;
                while MMloglik > 400000
                    [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, ssData, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                    itercount = itercount + 1;
                end
                disp('    id        setsize   SOA       Deviance  SD(deg)   P(guess)  P(trans)   P(CW)');
                disp([id, setsize, soa, MMloglik/1000, MMparms]);
                MMSD(id, soa, setsize) = MMparms(1);
                MMguessing(id, soa, setsize) = MMparms(2);
                if npar > 2, MMtranspos(id, soa, setsize) = MMparms(3); end
                if npar > 3, MMcwattraction(id, soa, setsize) = MMparms(4); end
                MMpm(id, soa, setsize) = 1-sum(MMparms(2:end));
            else
                disp('    id        setsize   SOA       mean(error)');
                disp([id, setsize, soa, Mdevobs(id, soa, setsize)]);
            end
            
        end % for cueing
    end %for setsize
    
    %fit IM (to all set size conditions jointly)
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

% Plot Mean(Deviation) as functions of set size

% Stand-Alone Figure of Precision (for direct comparison with Bays et al.)
PreFigure;
plotvector = squeeze(mean(Mprecision,1));  % precision in radians
plot(SOA, plotvector);
PostFigure([0, max(SOA), 0, 10], 'SOA', 'Precision (rad)', [], vec2legend(1:E.maxsetsize));

% Stand-Alone Figure of Errors
PreFigure;
plotvector = squeeze(mean(Mdevobs,1));
plot(SOA, plotvector);
PostFigure([0, max(SOA), 0, 45], 'SOA', 'Error (deg)', [], vec2legend(1:E.maxsetsize));



% Plot Mixture Model Parameters over Setsize
if fitMM
    meanMMpm = squeeze(mean(MMpm,1));
    meanK = bsxfun(@times, meanMMpm, 1:E.maxsetsize);
    PreFigure([], [], 2);
    subplot(3,2,1);
    plot(SOA, squeeze(mean(MMSD,1)));
    PostFigure([0, max(SOA), 0, max(max(mean(MMSD,1)))+0.5], 'SOA', 'Mean SD', 'SD from Bays Mixture', vec2legend(1:E.maxsetsize));
    subplot(3,2,2);
    plot(SOA, meanK');
    PostFigure([0, max(SOA), 0, 6], 'SOA', 'K', 'K from Bays Mixture', vec2legend(1:E.maxsetsize));
    subplot(3,2,3);
    plot(SOA, squeeze(mean(MMpm,1)));
    PostFigure([0, max(SOA), 0, 1], 'SOA', 'Mean P(m)', 'P(mem)');
    subplot(3,2,4);
    plot(SOA, squeeze(mean(MMguessing,1)));
    PostFigure([0, max(SOA), 0, 1], 'SOA', 'Mean P(guess)', 'P(guess)');
    subplot(3,2,5);
    plot(SOA, squeeze(mean(MMtranspos,1)));
    PostFigure([0, max(SOA), 0, 0.5], 'SOA', 'Mean P(trans)', 'P(transpos)');
    subplot(3,2,6);
    plot(SOA, squeeze(mean(MMcwattraction,1)));
    PostFigure([0, max(SOA), 0, 0.5], 'SOA', 'Mean P(wheel)', 'P(wheel attraction)');
end

if fitIM
    disp('      b       a          s       P.kappa      P.kappaf      r');
    disp(mean(IMparms, 1));
    disp('      b       a          s       P.kappa      P.kappaf      r');
    disp(std(IMparms, 1));
end

D.Mdevobs = Mdevobs;
if fitMM
    D.MMSD = MMSD;
    D.MMpm = MMpm;
    D.MMguessing = MMguessing;
    D.MMtranspos = MMtranspos;
    D.MMcwattraction = D.MMcwattraction;
end

%%% Save results
if E.saveResults == 1
    fid = fopen(['IMSim.MaskSOA.dat'], 'w');
    for id = 1:E.nsubj
        for setsize = 1:E.maxsetsize
            for soa = 1:length(SOA)
                fprintf(fid, '%d %d %d %d  ', id, setsize, soa, Mdevobs(id, soa, setsize));
                for ii = 1:length(C.indVar)
                    fprintf(fid, '%d ', ParX(id, ii));
                end                
                if (fitMM == 1), fprintf(fid, '%d %d %d %d', MMtranspos(id, soa, setsize), MMguessing(id, soa, setsize), MMcwattraction(id, soa, setsize), MMSD(id, soa, setsize)); end
                fprintf(fid, '\n');
            end
        end
    end
    fclose(fid);
end




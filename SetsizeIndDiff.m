% Simulation of Set-size and pre-cue/retro-cue effects

function [Kestimate, MMSD, MMguessing, MMtranspos, mCDAeeg, GroupPar] = SetsizeIndDiff(Model, fitIM)

global P
global E
global C

E.PreRetro = 2;
E.maxsetsize = 8;

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end


% Extract group means of parameters
if ~isempty(C.indVar)
    for ii = 1:length(C.indVar)
        eval(['meanIndVar(ii) = P.', C.indVar{ii}, ';']);  % use the parameter value in P as the mean 
    end
end
if ~isempty(C.groupVar)
    if ismember(C.groupVar, fields(P))
        groupVarIdx = find(strcmp(C.groupVar, C.indVar));
        eval(['parVal = P.', C.groupVar, ';']);
        meanGroupVar = parVal + [-2, -1, 0, 1, 2]*C.SDfactor(groupVarIdx)*parVal;
        sdGroupVar = C.SDfactor(groupVarIdx);
        GroupPar = zeros(E.ngroups, E.nsubj); 
    else
        meanGroupVar = 0;
        sdGroupVar = 0;
        groupVarIdx = 0; 
        GroupPar = 0; 
    end
    E.ngroups = length(meanGroupVar);
else
    E.ngroups = 1;
end

fitMM = 1;
option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

% initialize parameters of mixture model, IM, observed mean Deviation

Mdevobs = NaN(E.nsubj, E.ngroups, E.maxsetsize);
Mrt = NaN(E.nsubj, E.ngroups, E.maxsetsize);
Mgate = NaN(E.nsubj, E.ngroups, E.maxsetsize);
MMSD = NaN(E.nsubj, E.ngroups, E.maxsetsize);
MMguessing = NaN(E.nsubj, E.ngroups, E.maxsetsize);
MMtranspos = NaN(E.nsubj, E.ngroups, E.maxsetsize);
MMcwattraction = NaN(E.nsubj, E.ngroups, E.maxsetsize);
MMPmem = NaN(E.nsubj, E.ngroups, E.maxsetsize);
mCDAeeg = NaN(E.nsubj, E.ngroups, E.maxsetsize);


[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, E.maxsetsize), 1:360);  %Colors = E.ntrials x E.maxsetsize x [1:360]

P.nbmean = P.nb;  % keep the mean of nb for normalization of connection-weight matrices
Parms = zeros(E.ngroups*E.nsubj, length(C.indVar));

for group = 1:E.ngroups
    
    for ii = 1:length(C.indVar)
        if C.logistVar(ii)==1
            maxVar = logit(C.maxIndVar(ii)); meanVar = logit(meanIndVar(ii)); 
            logitPar = min(maxVar, randn(1, E.nsubj) + meanVar );  % SDfactor is ignored here because for parameters close to 0.5, meanVar gets close to 0, which would lead to zero variance. We just set SD=1 on the logit scale
            IndPar(ii,:) = logist(logitPar);
        else
            IndPar(ii,:) = min(C.maxIndVar(ii), max(0, randn(1, E.nsubj)*C.SDfactor(ii)*meanIndVar(ii) + meanIndVar(ii) ));
        end
        Parms(((group-1)*E.nsubj+1):group*E.nsubj, ii) = IndPar(ii,:);
    end
    if groupVarIdx > 0
        GroupPar(group,:) = min( C.maxGroupVar(groupVarIdx), max(0,  randn(1, E.nsubj)*sdGroupVar + meanGroupVar(group) ) );
        Parms(((group-1)*E.nsubj+1):group*E.nsubj, groupVarIdx) = GroupPar(group,:); % replace the individual parameter values for the grouping parameter by the group's values
    end
    
    plotx = ceil(sqrt(length(C.indVar)));
    ploty = ceil(sqrt(length(C.indVar)));
    PreFigure;
    for ii = 1:length(C.indVar)
        subplot(plotx, ploty, ii);
        histogram(IndPar(ii,:));
        PostFigure([], C.indVar{ii});
    end
    
    Array = zeros(E.nsubj*E.maxsetsize*E.ntrials, E.maxsetsize);
    Target = zeros(1, E.nsubj*E.maxsetsize*E.ntrials);
    Response = zeros(1, E.nsubj*E.maxsetsize*E.ntrials);
    tcount = 1; %trial count
    
    for id = 1:E.nsubj
        
        for ii = 1:length(C.indVar)
            eval(['P.', C.indVar{ii}, ' = IndPar(ii, id);']);
        end
        if groupVarIdx > 0, eval(['P.', C.groupVar, ' = GroupPar(group, id);']); end % replace ind par by group par for grouping parameter
        P.nb = round(P.nb);
        
        % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
        CreateStimuli;
        CreateMapping(E.calibrateAmp==2);
        
        for setsize = 1:E.maxsetsize
            fdistance = zeros(1,E.ntrials);
            rt = zeros(1,E.ntrials);
            GateClosed = zeros(1,E.ntrials);
            EEG_G = zeros(E.ntrials, 1);
            EEG_W = zeros(E.ntrials, 1); 
            Probedpos = zeros(E.ntrials,1);
            Pangle = zeros(E.ntrials,E.maxsetsize);
            Cangle = zeros(E.ntrials,E.maxsetsize+1);
            Targ = zeros(E.ntrials,1);
            Resp = zeros(E.ntrials,1);
            Setsize = zeros(E.ntrials,1);
            
            for trial = 1:E.ntrials
                
                output = Model(P, setsize, 1);
                
                Array(tcount,:) = [output.F(1:setsize), zeros(1, E.maxsetsize-setsize)];
                Target(tcount) = output.F(1);
                Response(tcount) = output.response;
                rt(trial) = output.rt;
                fdistance(trial) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
                EEG_W(trial, 1) = sum(abs(output.wx(:)));
                EEG_G(trial, 1) = sum(output.g); 
                GateClosed(trial) = mean(output.g); 
                tcount = tcount+1;
                
                %collect data for further modeling
                Probedpos(trial) = output.L(1);
                Pangle(trial,:) = output.L(1:E.maxsetsize);
                Cangle(trial,1:setsize) = output.F(1:setsize);
                Cangle(trial,E.maxsetsize+1) = output.CWcolor;
                Targ(trial) = output.F(1);
                Resp(trial) = output.response;
                Setsize(trial) = setsize;
                
            end
            
            Mdevobs(id, group, setsize) = mean(abs(fdistance));  %mean deviation
            Mrt(id, group, setsize) = mean(rt);
            Mgate(id, group, setsize) = mean(GateClosed);
            if C.CDA == 1, mCDAeeg(id, group, setsize) = mean(EEG_G); end
            if C.CDA == 2, mCDAeeg(id, group, setsize) = mean(EEG_W); end
            
            ssData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting
            ssData.cueing = repmat(1, length(Response), 1);
            ssData.preretro = E.PreRetro;
            ssD.setsize = Setsize;
            ssD.response = Resp;
            ssD.L = round(C.Location(Pangle));
            ssD.Color = Cangle;
            ssD.cueing = repmat(1, length(Response), 1);
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
                npar = 2;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture, 4 -> Souza & Oberauer mixture (iC.ncluding color-wheel attraction)
                MMloglik = 500000;
                itercount = 0;
                while MMloglik > 400000
                    [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, ssData, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                    itercount = itercount + 1;
                end
                varnames = {'id', 'setsize' 'iter', 'loglik/1000', 'SD', 'P(guess)', 'P(swap)', 'P(wheel)'};
                parXT = array2table(round([id, setsize, itercount, MMloglik/1000, MMparms], 2), 'VariableNames', varnames(1:(4+npar)));
                disp(parXT);
                MMSD(id, group, setsize) = MMparms(1);
                MMguessing(id, group, setsize) = MMparms(2);
                if npar > 2, MMtranspos(id, group, setsize) = MMparms(3); end
                if npar > 3, MMcwattraction(id, group, setsize) = MMparms(4); end
                MMPmem(id, group, setsize) = 1-MMparms(2:npar);
                
            end
            
        end %for setsize
        
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
            disp('     id       iter     loglik/1000  B       A       s        kappa       kappaf      cred');
            disp([id, itercount, IMloglik/1000, IMparms(id,:)]);
            pred = IM(IMparms(id,:), Data, 1);
            Dev = abs(wrap(repmat(1:360, size(Data.response,1), 1) - repmat(Data.response, 1, 360), 180));
            predDev = sum(Dev .* pred, 2); % Weights deviation from each possible angle with that angle's predicted probability (of being the response). Computes weighted average deviation for predicted SD
            predMDevIM(id, :) = aggregate(Data.setsize, predDev);
        end
        

    end  % for ID
    
    ResponseDistrib(Array, Target, Response);
    
end % for group

% Plot Mean(Deviation) as functions of set size

legendtext = vec2legend(meanGroupVar);
PreFigure;
subplot(2,2,1);
plotvector = squeeze(mean(Mdevobs,1))';
plot(plotvector);
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Deviation', 'Mean Deviation', legendtext);
subplot(2,2,2);
plotvector = squeeze(mean(Mgate,1))';
plot(plotvector);
PostFigure([0.8, E.maxsetsize+0.2, 0, 1], 'Setsize', 'Deviation', 'Prop. Binding Recruited', legendtext);
subplot(2,2,3);
plotvector = squeeze(mean((MMPmem),1))';
plot(plotvector);
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'P(mem)', 'P(mem) from Mixture Model', legendtext);
subplot(2,2,4);
plotvector = squeeze(mean(MMSD,1))';
plot(plotvector);
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'SD', 'SD from Mixture Model', legendtext);

% remove extreme outliers from mCDAeeg

mCDAeeg(mCDAeeg > 4*median(mCDAeeg(:))) = NaN; 

% put the groups together into a single sample
MCDA_EEG = [];
MGuess = [];
for group = 1:E.ngroups
    groupCDA = squeeze(mCDAeeg(:,group,:));
    groupGuessing = squeeze(MMguessing(:,group,:));
    groupParms = Parms(((group-1)*E.nsubj+1):(group*E.nsubj), :);
    [x,~] = find(isnan(groupCDA));
    if ~isempty(x)
        groupCDA(x,:) = [];
        groupGuessing(x,:) = [];
        groupParms(x,:) = [];
    end
    MCDA_EEG = [MCDA_EEG; groupCDA];
    MGuess = [MGuess; groupGuessing];
end

PreFigure;
plotvector = squeeze(nanmean(mCDAeeg,1))';
plot(plotvector, '-r', 'LineWidth', 5);
hold on
plot(MCDA_EEG', '-', 'color',[0,0,0]+0.5);
if C.CDA == 1, title = 'Absolute Binding-Pool Activation'; end
if C.CDA == 2, title = 'Gating Activation'; end
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*max(max(MCDA_EEG))], 'Setsize', 'CDA', title, legendtext);

eCDA0 = MCDA_EEG(:,4);  %weight-based -> EEG CDA at set size 4
eCDA14 = MCDA_EEG(:,4) - MCDA_EEG(:,1);
eCDA24 = MCDA_EEG(:,4) - MCDA_EEG(:,2);
Kestimate = zeros(size(MGuess,1), E.maxsetsize); 
for k = 1:E.maxsetsize
    Kestimate(:,k) = k * (1-MGuess(:,k));
end
K = mean(Kestimate(:, 4:E.maxsetsize), 2);


PreFigure;
subplot(2,2,1);
plot(K, eCDA0, 'o');
PostFigure([0, 6, 0, max(eCDA0)*1.05], 'K', 'CDA', '@Setsize 4');
subplot(2,2,2);
plot(K, eCDA14, 'o');
PostFigure([0, 6, min(eCDA14), max(eCDA14)*1.05], 'K', 'CDA', '@Setsize 4 - @Setsize 1');
subplot(2,2,3);
plot(K, eCDA24, 'o');
PostFigure([0, 6, min(eCDA24), max(eCDA24)*1.05], 'K', 'CDA', '@Setsize 4 - @Setsize 2');

% plot aCDA24 separately for PowerPoint presentations
PreFigure;
plot(K, eCDA24, 'o');
PostFigure([0, 5, min(eCDA24), max(eCDA24)*1.05], 'Memory Capacity (K)', 'CDA@4 - CDA@2');

if fitIM
    disp('          B        A       s       kappa       kappaf       Cred');
    disp(median(IMparms));
end

% Correlation matrix

ParmsPlus = [groupParms, K, eCDA0, eCDA14, eCDA24];
corrX = corrcoef(ParmsPlus);
varnames = [C.indVar, {'K', 'CDA4', 'CDA4_1', 'CDA4_2'}];
corrXT = array2table(round(corrX, 2), 'VariableNames', varnames, ...
    'RowNames', varnames);
disp(corrXT);


%%% Save results

save IMSim.SetsizeIndDiff.mat

if E.saveResults == 1
    fid = fopen('IMSim.SetsizeIndDiff.dat', 'w');
    for group = 1:E.ngroups
        for id = 1:E.nsubj
            for setsize = 1:E.maxsetsize
                fprintf(fid, '%d %d %d %d %d  ', id, setsize, group, Mdevobs(id, group, setsize), mCDAeeg(id, group, setsize));
                if (fitMM == 1), fprintf(fid, '%d %d %d', MMtranspos(id, group, setsize), MMguessing(id, group, setsize), MMSD(id, group, setsize)); end
                fprintf(fid, '\n');
            end
        end
    end
    fclose(fid);
    
    % save individual parameters
    fid = fopen('IMSim.SetsizeIndDiff.Parms.dat', 'w');
    for ii = 1:length(C.indVar), fprintf(fid, '%s ', C.indVar{ii}); end   % parameter names as column labels in first row
    fprintf(fid, '\n');
    for id = 1:size(Parms,1)  % goes across all groups and individuals within groups
        for ii = 1:length(C.indVar)
            fprintf(fid, '%d ', Parms(id, ii));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    
end



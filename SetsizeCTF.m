function [MMSD, MMguessing, MMtranspos] = SetsizeCTF(Model, fitMM, fitIM)
% Simulation of set-size effect on success of continuous-reproduction of locations, CDA,
% and reconstruction of channel tuning functions (CTF) through an inverted
% encoding model.

global P
global E
global C

items4CTF = 1;  % 1: all array items, sorted by input position, 2 = only target

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

% set up the IEM
CreateStimuli;
nElectrodes = 50;
eRegular = 0;  % factor of regular (distance-graded) to random projection weights in eW
eGrad = 1;     % generalization gradient in space for regular mapping of screen location (of the stimulus) to head location of electrode responding to it
eNoise = 0.25;  % trial-by-trial noise added to the EEG signal
nChannels = 9;  % an odd number has the advantage that there is a channel at the center of the scale, on which all CTFs are centered. 

E.PreRetro = 2;
E.maxsetsize = 6;
option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

% initialize parameters of mixture model, IM, observed mean Deviation

Mdevobs = NaN(E.nsubj, E.maxsetsize);
Mrt = NaN(E.nsubj, E.maxsetsize);
MMSD = NaN(E.nsubj, E.maxsetsize);
MMguessing = NaN(E.nsubj, E.maxsetsize);
MMtranspos = NaN(E.nsubj, E.maxsetsize);
MMcwattraction = NaN(E.nsubj, E.maxsetsize);
MMPmem = NaN(E.nsubj, E.maxsetsize);

[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, E.maxsetsize), 1:360);  %Colors = E.ntrials x E.maxsetsize x [1:360]
Array = zeros(E.nsubj*E.maxsetsize*E.ntrials, E.maxsetsize);
Target = zeros(1, E.nsubj*E.maxsetsize*E.ntrials);
Response = zeros(1, E.nsubj*E.maxsetsize*E.ntrials);
tcount = 1; %trial count

for id = 1:E.nsubj
    
    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end
    
    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);
    [basisSet, eW, eW2, nElectrodes, channelCenters] = CreateIEM(nElectrodes, nChannels, eRegular, eGrad);
    
    % Train IEM with set size = 1
    setsize = 1;
    nfactor = 1;
    EEG_WX = zeros(nfactor*E.ntrials, nElectrodes);
    EEG_FX = zeros(nfactor*E.ntrials, nElectrodes);
    StimMaskLoc = zeros(nfactor*E.ntrials, C.nc);
    StimMaskFeat = zeros(nfactor*E.ntrials, C.nstim);
    for trial = 1:(nfactor*E.ntrials)
        output = Model(P, setsize, 1);  % cueing = 1 (no cue)
        %EEG_WX(trial, :) = output.wx * eW2 + randn(1,nElectrodes)*eNoise; % read out w into electrodes directly
        EEG_WX(trial,:) = (((output.context * output.wx(1:C.nLocCat, :)) * output.wx((C.nLocCat+1):end, :)') * C.Mapping') * eW + randn(1,nElectrodes)*eNoise;  % feed last-used context into weight matrix -> reactivate content -> project onto electrodes
        EEG_FX(trial,:) = output.SpatAttn' * eW + randn(1,nElectrodes)*eNoise; % read out locations (averaging over features) from feature map
        StimMaskLoc(trial, round(C.Location(output.L(1:setsize)))) = 1; % stimulus mask: codes the stimulus location (set to 1 at presented location(s), and 0 everywhere else)
        StimMaskFeat(trial, round(C.feature(output.F(1:setsize)))) = 1; % stimulus mask: codes the stimulus feature (set to 1 at presented feature(s), and 0 everywhere else)
    end
    WwxLoc = TrainIEM(StimMaskLoc, basisSet, EEG_WX);    % train IEM
    WfxLoc = TrainIEM(StimMaskLoc, basisSet, EEG_FX);    % train IEM
    WwxFeat = TrainIEM(StimMaskFeat, basisSet, EEG_WX);    % train IEM
    
    for setsize = 1:E.maxsetsize
        
        fdistance = zeros(1,E.ntrials);
        rt = zeros(1,E.ntrials);
        Probedpos = zeros(E.ntrials,1);
        Pangle = zeros(E.ntrials,E.maxsetsize);
        Cangle = zeros(E.ntrials,E.maxsetsize+1);
        Targ = zeros(E.ntrials,1);
        Resp = zeros(E.ntrials,1);
        Setsize = zeros(E.ntrials,1);
        EEG_WX = zeros(E.ntrials, nElectrodes);
        EEG_FX = zeros(E.ntrials, nElectrodes);
        ItemIdx = zeros(E.ntrials, setsize);
        
        for trial = 1:E.ntrials
            
            output = Model(P, setsize, 1);
            
            % collect data for error-distribution plotting
            Array(tcount,:) = [output.F(1:setsize), zeros(1, E.maxsetsize-setsize)];
            Target(tcount) = output.F(1);
            Response(tcount) = output.response;
            tcount = tcount+1;
            
            %collect data for further modeling
            rt(trial) = output.rt;
            fdistance(trial) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
            Probedpos(trial) = output.L(1);
            Pangle(trial,:) = output.L(1:E.maxsetsize);
            Cangle(trial,1:setsize) = output.F(1:setsize);
            Cangle(trial,E.maxsetsize+1) = output.CWcolor;
            Targ(trial) = output.F(1);
            Resp(trial) = output.response;
            Setsize(trial) = setsize;
            %EEG_WX(trial, :) = output.wx * eW2 + randn(1,nElectrodes)*eNoise; % read out w into electrodes directly
            EEG_WX(trial,:) = (((output.context * output.wx(1:C.nLocCat, :)) * output.wx((C.nLocCat+1):end, :)') * C.Mapping') * eW + randn(1,nElectrodes)*eNoise;  % feed last-used context into weight matrix -> reactivate content -> project onto electrodes
            EEG_FX(trial,:) = output.SpatAttn' * eW + randn(1,nElectrodes)*eNoise; % read out locations (averaging over features) from feature map
            Inpos = output.Inpos; % returns the item index in their (random) order of encoding; item index 1 is the target!
            if items4CTF == 1
                for item = 1:setsize
                    ItemIdx(trial,Inpos(item)) = item; % sort items by their input position
                end
            end
        end
        
        if items4CTF == 2, ItemIdx = ones(E.ntrials,1); end   % ones: use only the target
        
        %CTF(id,setsize).meanCTF_WX = IEM(StimMask, basisSet, EEG_WX, Pangle, ItemIdx);
        %CTF(id,setsize).meanCTF_FX = IEM(StimMask, basisSet, EEG_FX, Pangle, ItemIdx);
        
        CTF(id,setsize).meanCTF_WX_Loc = ApplyIEM(WwxLoc, EEG_WX, 360*Pangle/C.nloc, channelCenters, ItemIdx);
        CTF(id,setsize).meanCTF_FX_Loc = ApplyIEM(WfxLoc, EEG_FX, 360*Pangle/C.nloc, channelCenters, ItemIdx);
        CTF(id,setsize).meanCTF_WX_Feat = ApplyIEM(WwxFeat, EEG_WX, Cangle, channelCenters, ItemIdx);
        
        Mdevobs(id, setsize) = mean(abs(fdistance));  %mean deviation
        Mrt(id, setsize) = mean(rt);
        
        ssData = Dataprocessing(Probedpos, Pangle, Cangle, Targ, Resp, Setsize, Colorgrid);   %prepare data for model fitting
        ssData.cueing = repmat(1, E.ntrials, 1);
        ssData.preretro = E.PreRetro;
        ssD.setsize = Setsize;
        ssD.response = Resp;
        ssD.L = round(C.Location(Pangle));
        ssD.Color = Cangle;
        ssD.cueing = repmat(1, E.ntrials, 1);
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
            disp([id, setsize, itercount, MMloglik/1000, MMparms]);
            MMSD(id, group, setsize) = MMparms(1);
            MMguessing(id, group, setsize) = MMparms(2);
            if npar > 2, MMtranspos(id, group, setsize) = MMparms(3); end
            if npar > 3, MMcwattraction(id, group, setsize) = MMparms(4); end
            MMPmem(id, group, setsize) = 1-MMparms(2:npar);
            
        end
        
        disp('    ID        Setsize   Error      ');
        disp([id, setsize, Mdevobs(id, setsize)]);
        
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

% plot CTFs for location from W

if E.saveResults == 1, fid = fopen('IMSim.SetsizeCTF.CTF.W2loc.dat', 'w'); end
PreFigure([], [], 2);
for setsize = 1:E.maxsetsize
    nItems = 1 + (2-items4CTF)*(setsize-1);
    mCTF = zeros(nItems, nChannels);
    for id = 1:E.nsubj
        for item = 1:nItems
            mCTF(item,:) = mCTF(item,:) + CTF(id,setsize).meanCTF_WX_Loc(item,:);
        end
    end
    mCTF = mCTF./E.nsubj;
    subplot(2,3,setsize);
    plot(channelCenters-180, mCTF);
    PostFigure([-180, 180, 0, 2], 'Location', 'CTF Response from W', ['Setsize: ', mat2str(setsize)], vec2legend(1:nItems));
end
if (E.saveResults == 1), fclose(fid); end

% plot CTFs for feature from W

if E.saveResults == 1, fid = fopen('IMSim.SetsizeCTF.CTF.W2feat.dat', 'w'); end
PreFigure([], [], 2);
for setsize = 1:E.maxsetsize
    nItems = 1 + (2-items4CTF)*(setsize-1);
    mCTF = zeros(nItems, nChannels);
    for id = 1:E.nsubj
        for item = 1:nItems
            mCTF(item,:) = mCTF(item,:) + CTF(id,setsize).meanCTF_WX_Feat(item,:);
        end
    end
    mCTF = mCTF./E.nsubj;
    subplot(2,3,setsize);
    plot(channelCenters-180, mCTF);
    PostFigure([-180, 180, 0, 2], 'Feature Value', 'CTF Response from W', ['Setsize: ', mat2str(setsize)], vec2legend(1:nItems));
end
if (E.saveResults == 1), fclose(fid); end

% plot CTFs for locations from feature maps

if E.saveResults == 1, fid = fopen('IMSim.SetsizeCTF.CTF.FX2loc.dat', 'w'); end
PreFigure([], [], 2);
for setsize = 1:E.maxsetsize
    nItems = 1 + (2-items4CTF)*(setsize-1);
    mCTF = zeros(nItems, nChannels);
    for id = 1:E.nsubj
        for item = 1:nItems
            mCTF(item,:) = mCTF(item,:) + CTF(id,setsize).meanCTF_FX_Loc(item,:);
        end
    end
    mCTF = mCTF./E.nsubj;
    subplot(2,3,setsize);
    plot(channelCenters-180, mCTF);
    PostFigure([-180, 180, 0, 2], 'Location', 'CTF Response from FX', ['Setsize: ', mat2str(setsize)], vec2legend(1:nItems));
    if E.saveResults == 1
        for item = 1:nItems
            fprintf(fid, '%d %d  ', setsize, item);
            for channel = 1:nChannels
                fprintf(fid, '%d ', mCTF(item, channel));
            end
            fprintf(fid, '\n');
        end
    end
end
if (E.saveResults == 1), fclose(fid); end


%%% Save results
if E.saveResults == 1
    fid = fopen('IMSim.SetsizeCTF.Behav.dat', 'w');
    for id = 1:E.nsubj
        for setsize = 1:E.maxsetsize
            fprintf(fid, '%d %d   %d %d  ', id, setsize, Mdevobs(id, setsize), Mrt(id, setsize));
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
end

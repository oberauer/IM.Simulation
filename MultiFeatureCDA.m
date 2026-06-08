function D = MultiFeatureCDA(Model, maxSetsize)
% Simulation of Set-size and number of features, on Alpha suppression and CDA, with continuous reproduction

global P
global E
global C

E.maxsetsize = maxSetsize;
E.presentation = 2;  % sequential presentation (just for fun)

maxFeat = 4;    % maximum number of features (or feature dimensions) per object
items4CTF = 1;  % 1: all array items, sorted by input position, 2 = only target

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

% set up the IEM
C.nfeatures = maxFeat;   % for creating stimuli
CreateStimuli;
nElectrodes = 50;
eRegular = 0;  % factor of regular (distance-graded) to random projection weights in eW
eGrad = 1;     % generalization gradient in space for regular mapping of screen location (of the stimulus) to head location of electrode responding to it
eNoise = 0.25;  % trial-by-trial noise added to the EEG signal
nChannels = 9;

Mdevobs = NaN(E.nsubj, C.nfeatures, E.maxsetsize, E.maxsetsize);  % id, nfeat, inpos, setsize
mCDAg = NaN(E.nsubj, C.nfeatures, E.maxsetsize);  % id, nfeat, setsize
mCDAw = NaN(E.nsubj, C.nfeatures, E.maxsetsize);  % id, nfeat, setsize
mAlpha = NaN(E.nsubj, C.nfeatures, E.maxsetsize);  % id, nfeat, setsize
Amplify = NaN(E.nsubj, 1);

[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, E.maxsetsize), 1:360);  %Colors = E.ntrials x E.maxsetsize x [1:360]

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
    Amplify(id) = C.amplify;
    
    for nfeat = 1:maxFeat
        
        E.nfeat = nfeat;
        C.nfeatures = nfeat;
        
        for ssidx = 0:E.maxsetsize
            
            setsize = max(1, ssidx);
            fdistance = NaN(E.ntrials, setsize);  % distance (target, response) for each input position (position in presentation order)
            Pangle = zeros(E.ntrials,E.maxsetsize);
            fx = zeros(1,E.ntrials);
            EEG_G = zeros(E.ntrials, 1); 
            EEG_W = zeros(E.ntrials, 1); 
            EEG_FX = zeros(E.ntrials, nElectrodes);
            StimMask = zeros(E.ntrials, C.nc);
            ItemIdx = zeros(E.ntrials, setsize);
            
            for trial = 1:E.ntrials
                
                output = Model(P, setsize, 1);  % cueing = 1 (no cue)
                fdistance(trial, output.Inpos(1)) = wrap(output.response-output.F(1), 180);   %calculate distance between response and true feature in feature space (degrees!)
                Pangle(trial,:) = output.L(1:E.maxsetsize);
                fx(trial) = sum(sum(output.fx));   % sum of activation in feature map of target feature dimension -> Alpha power
                EEG_W(trial) = sum(abs(output.wx(:)));
                EEG_G(trial) = sum(output.g); % closed gating units -> CDA
                EEG_FX(trial, :) = output.SpatAttn' * eW + randn(1,nElectrodes)*eNoise; % attended locations
                StimMask(trial, round(C.Location(output.L(1:setsize)))) = 1; % stimulus mask: codes the stimulus location (set to 1 at presented location(s), and 0 everywhere else)
                Inpos = output.Inpos; % returns the item index in their (random) order of encoding; item index 1 is the target!
                if items4CTF == 1
                    for item = 1:setsize
                        ItemIdx(trial,Inpos(item)) = item; % sort items by their input position
                    end
                end
                tcount = tcount+1;
            end
            
            if items4CTF == 2, ItemIdx = ones(E.ntrials,1); end   % ones: use only the target
            
            if ssidx == 0, W = TrainIEM(StimMask, basisSet, EEG_FX); end   % train on set size 1 only
            if ssidx > 0, CTF(id,nfeat,setsize).meanCTF_FX = ApplyIEM(W, EEG_FX, 360*Pangle/C.nloc, channelCenters, ItemIdx); end
            
            for inpos = 1:setsize
                Mdevobs(id, nfeat, setsize, inpos) = nanmean(abs(fdistance(:,inpos)));  %mean deviation 
                %MeanPreDeviation(id, nfeat, setsize, inpos) = mean(PreDeviation(:, inpos)); 
            end
            mCDAg(id, nfeat, setsize) = mean(EEG_G); 
            mCDAw(id, nfeat, setsize) = mean(EEG_W); 
            mAlpha(id, nfeat, setsize) = mean(sum(EEG_FX,2));
            
            disp('    ID        N(feat)   Setsize   Error     CDA/1000  Alpha      ');
            disp([id, nfeat, setsize, nanmean(Mdevobs(id, nfeat, setsize, :)), mCDAg(id, nfeat, setsize)/1000, mAlpha(id, nfeat, setsize)]);
            
        end %for setsize
        
    end % for nfeat
    
    
end  % for ID


% Plot Mean(Deviation) as function of set size and input position, with a
% separate plot for each number of features

PreFigure;
legendtext = {'In=1', 'In=2', 'In=3','In=4', 'In=5', 'In=6', 'In=7', 'In=8'};
for nfeat = 1:maxFeat
    subplot(2,2,nfeat);
    plotvector = squeeze(nanmean(Mdevobs(:,nfeat,:,:),1));  % select nfeat presentation, average over subjects (1)
    plot(plotvector);
    PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Deviation (Deg)', ['N(features) = ', mat2str(nfeat)], legendtext);
end
PreFigure;
legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
for nfeat = 1:maxFeat
    plotvector = squeeze(nanmean(Mdevobs(:,nfeat,:,:),1));  % select nfeat, average over subjects (1)
    subplot(2,2,nfeat);
    plot(plotvector');
    PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Input Position', 'Deviation (Deg)', ['N(features) = ', mat2str(nfeat)], legendtext);
end

% plot Mean(Deviation) as function of set size and number of features
PreFigure;
legendtext = {'NF=1', 'NF=2', 'NF=3', 'NF=4'};
plotvector = squeeze(nanmean(nanmean(Mdevobs, 4), 1)); % average over input position (4) and subjects (1)
plot(1:setsize, plotvector');
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Deviation (Deg)', [], legendtext);

% Plot CDA and Alpha Power as a function of set size and number of features

meanCDAg = zeros(4, E.maxsetsize);
meanCDAw = zeros(4, E.maxsetsize);
meanAlpha = zeros(4, E.maxsetsize);
for nfeat = 1:maxFeat
    meanCDAg(nfeat, :) = squeeze(mean(mCDAg(:,nfeat,:)));
    meanCDAw(nfeat, :) = squeeze(mean(mCDAw(:,nfeat,:)));
    meanAlpha(nfeat, :) = squeeze(mean(mAlpha(:,nfeat,:)));
end

PreFigure;
ymax = max(max(meanCDAg(:)), max(meanCDAw(:)));
subplot(1,2,1);
plot(1:E.maxsetsize, meanCDAg');
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*ymax], 'Set Size', 'CDA from Closed Gating Units', [], vec2legend(1:4));
subplot(1,2,2);
plot(1:E.maxsetsize, meanCDAg');
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*ymax], 'Set Size', 'CDA from abs(Weight Matrix)', [], vec2legend(1:4));

PreFigure;
ymax = max(max(meanAlpha));
plot(1:E.maxsetsize, meanAlpha');
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*ymax], 'Set Size', 'Alpha', [], vec2legend(1:4));


% plot CTFs from feature maps

SimSeq = {'1 Feature', '2 Features', '3 Features', '4 Features'};
for nfeat = 1:maxFeat
    PreFigure([], [], 2);
    for setsize = 1:E.maxsetsize
        nItems = 1 + (2-items4CTF)*(setsize-1);
        mCTF = zeros(nItems, nChannels);
        for id = 1:E.nsubj
            for item = 1:nItems
                mCTF(item,:) = mCTF(item,:) + CTF(id,nfeat,setsize).meanCTF_FX(item,:);
            end
        end
        mCTF = mCTF./E.nsubj;
        subplot(2,3,setsize);
        plot(channelCenters-180, mCTF);
        PostFigure([-180, 180, 0, 2], 'Feature Dimension', 'CTF of Location from FX', ['SS ', mat2str(setsize) ', ', SimSeq{nfeat}], vec2legend(1:nItems));
    end
end

D.Mdevobs = Mdevobs;
D.mCTF = mCTF;
D.CTF = CTF;
D.meanCDAg = meanCDAg;
D.meanCDAg = meanCDAw;
D.meanAlpha = meanAlpha;


%%% Save results

if E.saveResults == 1
    
    fid = fopen('IMSim.MultiFeatError.dat', 'w');
    for id = 1:E.nsubj
        for setsize = 1:E.maxsetsize
            for nfeat = 1:maxFeat
                for inpos = 1:setsize
                    fprintf(fid, '%d %d %d %d %d  ', id, setsize, nfeat, inpos, Mdevobs(id, nfeat, setsize, inpos));
                    fprintf(fid, '\n');
                end
            end
        end
    end
    fclose(fid);
    
    fid = fopen('IMSim.;ultiFeatCDA.dat', 'w');
    for id = 1:E.nsubj
        for setsize = 1:E.maxsetsize
            for nfeat = 1:4
                fprintf(fid, '%d %d %d %d %d  ', id, setsize, nfeat, mCDAg(id, nfeat, setsize), mAlpha(id, nfeat, setsize));
                fprintf(fid, '\n');
            end
        end
    end
    fclose(fid);
    
    fid = fopen('IMSim.MultiFeatCTF.dat', 'w');
    for id = 1:E.nsubj
        for setsize = 1:E.maxsetsize
            for nfeat = 1:maxFeat
                for item = 1:setsize
                    fprintf(fid, '%d %d %d %d    ', id, setsize, nfeat, item);
                    for channel = 1:nChannels
                        fprintf(fid, '%d ', CTF(id, nfeat, setsize).meanCTF_FX(item, channel));
                    end
                    fprintf(fid, '\n');
                end
            end
        end
    end
    fclose(fid);
    
end


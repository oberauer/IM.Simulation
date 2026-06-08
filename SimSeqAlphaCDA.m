function D = SimSeqAlphaCDA(Model, maxSetsize)
% Simulation of Set-size and simultaneous vs sequential presentation, on Alpha suppression and CDA, with continuous reproduction

global P
global E
global C

E.maxsetsize = maxSetsize;
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
nChannels = 9; 

Mdevobs = NaN(E.nsubj, 2, E.maxsetsize, E.maxsetsize, E.maxsetsize);  % id, simseq, setsize, inpos, outpos
Mstrength = NaN(E.nsubj, 2, E.maxsetsize, E.maxsetsize, E.maxsetsize);  % id, simseq, setsize, inpos, outpos
BindStrength = NaN(E.nsubj, E.maxsetsize, E.maxsetsize); % id, inpos, setsize; for sequential condition only (so far)
mCDA = NaN(E.nsubj, 2, E.maxsetsize);  % id, simseq, setsize
mAlpha = NaN(E.nsubj, 2, E.maxsetsize);  % id, simseq, setsize
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
    
    for simseq = 1:2
        
        E.presentation = simseq;
        if simseq == 1, E.prestime = 1; else, E.prestime = 0.5; E.ISI = 0; end
        
        for ssidx = 0:E.maxsetsize
            
            setsize = max(1, ssidx); 
            E.outsize = setsize;
            fdistance = zeros(E.ntrials, setsize, setsize);  % distance (target, response) for each input position (position in presentation order) x output position (position in test order)
            Pangle = zeros(E.ntrials,E.maxsetsize);
            fx = zeros(1,E.ntrials);
            EEG_W = zeros(E.ntrials, 1);
            EEG_G = zeros(E.ntrials, 1); 
            EEG_FX = zeros(E.ntrials, nElectrodes);
            StimMask = zeros(E.ntrials, C.nc); 
            ItemIdx = zeros(E.ntrials, setsize);
            encStrength = zeros(E.ntrials, setsize);
            bindStrength = zeros(E.ntrials, setsize);
            strengthIO = zeros(E.ntrials, setsize, setsize);   % encoding strength for each input position (position in presentation order) x output position (position in test order)
            
            for trial = 1:E.ntrials
                
                output = Model(P, setsize, 1);  % cueing = 1 (no cue)
                for outpos = 1:E.outsize
                    fdistance(trial, output.Inpos(outpos), outpos) = wrap(output.response(outpos)-output.F(outpos), 180);   %calculate distance between response and true feature in feature space (degrees!)
                    strengthIO(trial, output.Inpos(outpos), outpos) = output.Bstrength(output.Inpos(outpos)); 
                end
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
                encStrength(trial, :) = output.Strength;
                bindStrength(trial, :) = output.Bstrength';
                tcount = tcount+1;
            end
            
            if items4CTF == 2, ItemIdx = ones(E.ntrials,1); end   % ones: use only the target
        
            %CTF(id,simseq,setsize).meanCTF_WX = IEM(StimMask, basisSet, EEG_WX, Pangle, ItemIdx);  
            if ssidx == 0, W = TrainIEM(StimMask, basisSet, EEG_FX); end   % train on set size 1 only
            if ssidx > 0, CTF(id,simseq,setsize).meanCTF_FX = ApplyIEM(W, EEG_FX, 360*Pangle/C.nloc, channelCenters, ItemIdx); end
            
            for inpos = 1:setsize
                for outpos = 1:setsize
                    Mdevobs(id, simseq, setsize, inpos, outpos) = mean(abs(fdistance(:,inpos,outpos)))*setsize;  %mean deviation - need to multiply by setsize because matrix is setsize*setsize but has only setsize non-zero entries
                    Mstrength(id, simseq, setsize, inpos, outpos) = mean(strengthIO(:,inpos,outpos)); 
                end
            end
            if simseq == 2, BindStrength(id, setsize, 1:setsize) = mean(bindStrength); end
            if C.CDA == 1, mCDA(id, simseq, setsize) = mean(EEG_G); end
            if C.CDA == 2, mCDA(id, simseq, setsize) = mean(EEG_W); end
            mAlpha(id, simseq, setsize) = mean(sum(EEG_FX,2)); 
            
            disp('    ID        Simseq    Setsize   CDA/1000  Alpha      ');
            disp([id, simseq, setsize, mCDA(id, simseq, setsize)/1000, mAlpha(id, simseq, setsize)]);
            
        end %for setsize
        
    end % for simseq

    
end  % for ID

% Plot Mean(Deviation) as a function of set size and output position (simultaneous presentation)
PreFigure;
legendtext = {'Out=1', 'Out=2', 'Out=3','Out=4', 'Out=5', 'Out=6', 'Out=7', 'Out=8'};
plotvector = squeeze(nanmean(nanmean(Mdevobs(:,1,:,:,:),4),1));  % select simultaneous presentation, average over inpos (4) and subjects (1)
subplot(1,2,1);
plot(plotvector);
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Deviation (Deg)', 'Simultaneous', legendtext(1:setsize));
legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
subplot(1,2,2);
plot(plotvector');
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Output Position', 'Deviation (Deg)', 'Simultaneous', legendtext(1:setsize));


% Plot Mean(Deviation) as function of set size and output position (sequential presentation)
PreFigure;
legendtext = {'Out=1', 'Out=2', 'Out=3','Out=4', 'Out=5', 'Out=6', 'Out=7', 'Out=8'};
subplot(1,2,1);
plotvector = squeeze(nanmean(nanmean(Mdevobs(:,2,:,:,:),4),1));  % select sequential presentatino, average over inpos (4) and subjects (1)
plot(plotvector);
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Deviation (Deg)', 'Sequential', legendtext(1:setsize));

legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
subplot(1,2,2);
plot(plotvector');
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Output Position', 'Deviation (Deg)', 'Sequential', legendtext(1:setsize));


% Plot Mean(Deviation) as function of set size and input position (sequential presentation)
legendtext = {'In=1', 'In=2', 'In=3','In=4', 'In=5', 'In=6', 'In=7', 'In=8'};
PreFigure;
subplot(1,2,1);
plotvector = squeeze(nanmean(nanmean(Mdevobs(:,2,:,:,:),5),1));  % select sequential presentation, average over outpos (5) and subjects (1)
plot(plotvector);
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Deviation (Deg)', 'Sequential', legendtext(1:setsize));

legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
subplot(1,2,2);
plot(plotvector');
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Input Position', 'Deviation (Deg)', 'Sequential', legendtext(1:setsize));

% plot Mean(Deviation) as function of input x output position (sequential presentation, set size 6)

PreFigure;
legendtext = {'Out=1', 'Out=2', 'Out=3','Out=4', 'Out=5', 'Out=6', 'Out=7', 'Out=8'};
subplot(1,2,1);
plotvector = squeeze(nanmean(Mdevobs(:,2,6,:,:),1));  % select sequential presentation, set size 6, average over subjects (1)
plot(plotvector);
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Input Position', 'Deviation (Deg)', 'Sequential, SS=6', legendtext(1:setsize));
legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
subplot(1,2,2);
plot(squeeze(mean(BindStrength))'); 
PostFigure([0.8,setsize+0.2, 0, 1], 'Input Position', 'Binding Strength', 'Sequential', legendtext(1:setsize));

% Plot CDA and Alpha Power as a function of set size and simultaneous/sequential presentation 

CDAsim = squeeze(mean(mCDA(:,1,:), 1));
CDAseq = squeeze(mean(mCDA(:,2,:), 1));
Alphasim = squeeze(mean(mAlpha(:,1,:), 1));
Alphaseq = squeeze(mean(mAlpha(:,2,:), 1));

PreFigure;
ymax = max(max(CDAsim), max(CDAseq)); 
subplot(1,2,1);
plot(1:E.maxsetsize, CDAsim);  % simultaneous presentation 
if C.CDA == 1, CDAlegend = "CDA from closed gating units"; end
if C.CDA == 2, CDAlegend = "CDA from weight matrix"; end
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*ymax], 'Set Size', CDAlegend, 'Simultaneous');
subplot(1,2,2);
plot(1:E.maxsetsize, CDAseq);  % sequential presentation 
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*ymax], 'Set Size', CDAlegend, 'Sequential');

PreFigure;
ymax = max(max(Alphasim), max(Alphaseq)); 
subplot(1,2,1);
plot(1:E.maxsetsize, Alphasim);  % simultaneous presentation 
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*ymax], 'Set Size', 'Alpha', 'Simultaneous');
subplot(1,2,2);
plot(1:E.maxsetsize, Alphaseq);  % sequential presentation 
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*ymax], 'Set Size', 'Alpha', 'Sequential');


% plot CTFs from feature maps

SimSeq = {'Sim.', 'Seq.'};
for simseq = 1:2
    PreFigure([], [], 2);
    for setsize = 1:E.maxsetsize
        nItems = 1 + (2-items4CTF)*(setsize-1);
        mCTF = zeros(nItems, nChannels);
        for id = 1:E.nsubj
            for item = 1:nItems
                mCTF(item,:) = mCTF(item,:) + CTF(id,simseq,setsize).meanCTF_FX(item,:);
            end
        end
        mCTF = mCTF./E.nsubj;
        subplot(2,3,setsize);
        plot(channelCenters-180, mCTF);
        PostFigure([-180, 180, 0, 2], 'Feature Value', 'CTF Response from FX', ['SS ', mat2str(setsize) ', ', SimSeq{simseq}], vec2legend(1:nItems));
    end
end


% Plot MStrength as function of set size and output position (sequential presentation)
PreFigure;
legendtext = {'Out=1', 'Out=2', 'Out=3','Out=4', 'Out=5', 'Out=6', 'Out=7', 'Out=8'};
subplot(1,2,1);
plotvector = squeeze(nanmean(nanmean(Mstrength(:,2,:,:,:),4),1));  % select sequential presentation, average over inpos (4) and subjects (1)
plot(plotvector);
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Strength', 'Sequential', legendtext(1:E.maxsetsize));

legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
subplot(1,2,2);
plot(plotvector');
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Output Position', 'Strength', 'Sequential', legendtext(1:E.maxsetsize));

D.Mdevobs = Mdevobs;
D.mCDA = mCDA;
D.mAlpha = mAlpha;
D.CTF = CTF;

%%% Save results

if E.saveResults == 1
    
    fid = fopen('IMSim.SimSeqError.dat', 'w');
    for id = 1:E.nsubj
        for setsize = 1:E.maxsetsize
            for simseq = 1:2
                for inpos = 1:setsize
                    for outpos = 1:setsize
                        fprintf(fid, '%d %d %d %d %d %d  ', id, setsize, simseq, inpos, outpos, Mdevobs(id, simseq, setsize, inpos, outpos));
                        fprintf(fid, '\n');
                    end
                end
            end
        end
    end
    fclose(fid);
    
    fid = fopen('IMSim.SimSeqAlphaCDA.dat', 'w');
    for id = 1:E.nsubj
        for setsize = 1:E.maxsetsize
            for simseq = 1:2
                fprintf(fid, '%d %d %d %d %d  ', id, setsize, simseq, mCDA(id, simseq, setsize), mAlpha(id, simseq, setsize));
                fprintf(fid, '\n');
            end
        end
    end
    fclose(fid);
    
    fid = fopen('IMSim.SimSeqCTF.dat', 'w');
    for id = 1:E.nsubj
        for setsize = 1:E.maxsetsize
            for simseq = 1:2
                for item = 1:setsize
                    fprintf(fid, '%d %d %d %d    ', id, setsize, simseq, item);
                    for channel = 1:nChannels
                        fprintf(fid, '%d ', CTF(id, simseq, setsize).meanCTF_FX(item, channel));
                    end
                    fprintf(fid, '\n');
                end
            end
        end
    end
    fclose(fid);
    
end


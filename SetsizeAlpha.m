function [] = SetsizeAlpha(Model, maxSetsize)
% Simulation of Set-size, measuring time course of spatial attention
% modulation, extracting alpha power from it at each electrode as input to
% IEM

global P
global E
global C

C.nloc = 360;  % needs 360 locations
C.tstep = 0.005; % needs finer temporal resolution than usual
E.maxsetsize = maxSetsize;
E.presentation = 1;
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
Fs = 1/C.tstep;

% MdevByOut = NaN(E.nsubj, E.maxsetsize, E.maxsetsize);  % id, setsize, outpos
% MdevByIn = NaN(E.nsubj, E.maxsetsize, E.maxsetsize);  % id, setsize, inpos
Mdevobs = NaN(E.nsubj, E.maxsetsize, E.maxsetsize, E.maxsetsize);  % id, setsize, inpos, outpos
mAlpha = NaN(E.nsubj, E.maxsetsize);  % id, setsize
BindingStrength = NaN(E.nsubj, E.maxsetsize, E.maxsetsize); % id, setsize, serial position
NumberBound = NaN(E.nsubj, E.maxsetsize); % id, setsize


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
    
    for ssidx = 0:E.maxsetsize
        
        setsize = max(1, ssidx);
        E.outsize = setsize;
        fdistance = zeros(E.ntrials, setsize, setsize);  % distance (target, response) for each input position (position in self-chosen encoding order) x output position (position in test order)
        Pangle = zeros(E.ntrials,E.maxsetsize);
        EEG_FX = zeros(E.ntrials, nElectrodes);
        nBound = zeros(E.ntrials, 1); 
        bStrength = zeros(E.ntrials, setsize); 
        StimMask = zeros(E.ntrials, C.nc);
        ItemIdx = zeros(E.ntrials, setsize);
        
        for trial = 1:E.ntrials
            
            plotAlpha = id==1 && trial==1; 
            output = Model(P, setsize, plotAlpha);
            for outpos = 1:setsize
                %fdistance(trial, outpos) = wrap(output.response(outpos)-output.F(outpos), 180);   %calculate distance between response and true feature in feature space (degrees!)
                fdistance(trial, output.Inpos(outpos), outpos) = wrap(output.response(outpos)-output.F(outpos), 180);   %calculate distance between response and true feature in feature space (degrees!)
            end
            Pangle(trial,:) = output.L(1:E.maxsetsize);
            
            % Extract alpha power from the time course of attention
            sLength = size(output.SpatAttn,1); % signal length (number of time steps) = length of the FFT output
            for tIdx = 1:sLength
                Signal(tIdx,:) = output.SpatAttn(tIdx,:) * eW + randn(1,nElectrodes)*eNoise;
            end
            
            % Fourier transform
            Y = fft(Signal);
            Power2sided = abs(Y/sLength);
            Power = Power2sided(1:sLength/2+1,:);
            Power(2:end-1, :) = 2*Power(2:end-1, :);
            Hz = Fs*(0:(sLength/2))/sLength;
            AlphaFilteredSignal = eegfilt(Signal', Fs, 8, 12, sLength, 60)';
            AlphaInstantPower = zeros(size(Signal));
            for sensor = 1:nElectrodes
                AlphaInstantPower(:, sensor) = abs(hilbert(AlphaFilteredSignal(:,sensor)')).^2;
            end
            EEG_FX(trial, :) = mean(AlphaInstantPower(20:(sLength-20), :));  % average over all time points except a brief buffer at the beginning and the end, where instantaneous power is distorted
            
            if plotAlpha
                PreFigure([], [], 2);
                subplot(2,2,1);
                plot(Hz, mean(Power, 2));
                PostFigure([0, 1.05*max(Hz), 0, 1.1*max(mean(Power,2))], 'Hz', 'Power');
                subplot(2,2,2);
                plot(C.tstep:C.tstep:E.RI, AlphaFilteredSignal);
                PostFigure([0, E.RI, 0, 1.1*max(max(AlphaFilteredSignal))], 'Time (s)', 'Alpha Signal');
                subplot(2,2,3);
                plot(C.tstep:C.tstep:E.RI, mean(AlphaInstantPower, 2));
                PostFigure([0, E.RI, 0, 1.1*max(mean(AlphaInstantPower, 2))], 'Time (s)', 'Alpha Instantaneous Power');
            end
            
            StimMask(trial, round(C.Location(output.L(1:setsize)))) = 1; % stimulus mask: codes the stimulus location (set to 1 at presented location(s), and 0 everywhere else)
            Inpos = output.Inpos; % returns the item index in their (random) order of encoding; item index 1 is the target!
            if items4CTF == 1
                for item = 1:setsize
                    ItemIdx(trial,Inpos(item)) = item; % sort items by their input position
                end
            end
            bStrength(trial, :) = output.Bstrength(ItemIdx(trial,:))';  % sort the binding strengths by the item's input position
            nBound(trial) = length(nonzeros(unique(output.Bstrength)));

            tcount = tcount+1;
            
        end
        
        if items4CTF == 2, ItemIdx = ones(E.ntrials,1); end   % ones: use only the target
        
        if ssidx == 0, W = TrainIEM(StimMask, basisSet, EEG_FX); end   % train on set size 1 only
        if ssidx > 0, CTF(id,setsize).meanCTF_FX = ApplyIEM(W, EEG_FX, C.Location(Pangle), channelCenters, ItemIdx); end
        
%         for inpos = 1:setsize
%             MdevByIn(id, setsize, inpos) = mean(abs(fdistance(:, Inpos(inpos))));
%         end
%         for outpos = 1:setsize
%             MdevByOut(id, setsize, outpos) = mean(abs(fdistance(:,outpos)));
%         end
          for inpos = 1:setsize
             for outpos = 1:setsize
                 Mdevobs(id, setsize, inpos, outpos) = mean(abs(fdistance(:,inpos,outpos)))*setsize;  %mean deviation - need to multiply by setsize because matrix is setsize*setsize but has only setsize non-zero entries
             end
         end
        
        mAlpha(id, setsize) = mean(sum(abs(EEG_FX),2));
        BindingStrength(id, setsize, 1:setsize) = mean(bStrength);  
        NumberBound(id, setsize) = mean(nBound); 
        
        disp('    ID        Setsize   Alpha      ');
        disp([id, setsize, mAlpha(id, setsize)]);
        
    end %for setsize
    
    
end  % for ID

% Plot Mean(Deviation) as a function of set size and output position (simultaneous presentation)
PreFigure;
legendtext = {'Out=1', 'Out=2', 'Out=3','Out=4', 'Out=5', 'Out=6', 'Out=7', 'Out=8'};
plotvector = squeeze(nanmean(nanmean(Mdevobs, 3), 1));  % average over inpos (3) and subjects (1)
subplot(1,3,1);
plot(plotvector);
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*max(max(plotvector))], 'Setsize', 'Deviation (Deg)', [], legendtext);
legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
subplot(1,3,2);
plot(plotvector');
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*max(max(plotvector))], 'Output Position', 'Deviation (Deg)', [], legendtext);
subplot(1,3,3);
plotvector = squeeze(nanmean(nanmean(Mdevobs, 4), 1));  % average over outpos(4) and subjects (1)
plot(plotvector');
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*max(max(plotvector))], 'Input Position', 'Deviation (Deg)', [], legendtext);

PreFigure;
plot(1:E.maxsetsize, nanmean(mAlpha));
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*max(nanmean(mAlpha))], 'Set Size', 'Alpha Power');

% plot number of items bound at all, and binding strength
PreFigure;
legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
subplot(1,2,1);
plot(1:E.maxsetsize, mean(NumberBound));
PostFigure([0.8, E.maxsetsize+0.2, 0, E.maxsetsize], 'Setsize', 'N(Bound)', 'Number of Items Bound');
subplot(1,2,2);
plot(1:E.maxsetsize, squeeze(mean(BindingStrength))');
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.1*max(max(mean(BindingStrength)))], 'Input Position', 'Strength', 'Binding Strength', legendtext);

% plot CTFs from feature maps
PreFigure([], [], 2);
for setsize = 1:E.maxsetsize
    nItems = 1 + (2-items4CTF)*(setsize-1);
    mCTF = zeros(nItems, nChannels);
    for id = 1:E.nsubj
        for item = 1:nItems
            mCTF(item,:) = mCTF(item,:) + CTF(id,setsize).meanCTF_FX(item,:);
        end
    end
    mCTF = mCTF./E.nsubj;
    subplot(3,3,setsize);
    plot(channelCenters-180, mCTF);
    PostFigure([-180, 180, 0, 0.8], 'Feature Value', 'CTF Response from FX', ['SS ', mat2str(setsize)], vec2legend(1:nItems));
end


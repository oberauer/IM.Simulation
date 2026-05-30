function [] = TwoArrayISI
% Simulation of simultaneous presentation of two arrays with variable ISI
% (Jacob's experiment on consolidation)

global P
global E
global C

E.mask = 1;
E.prestime = 0.15;
E.MaskSOA = 0.15;  %
E.presentation = 1; % each array is presented simultaneously
E.outsize = 1;
InterArrayInterval = [0.2, 0.3, 0.4, 0.6, 1];

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

Mdevobs = NaN(E.nsubj, 3, 3, 2, length(InterArrayInterval));  % id, SS1, SS2, testedArray, III
Mbstrength = NaN(E.nsubj, 3, 3, 2, length(InterArrayInterval));  % id, SS1, SS2, Array, III

%[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x setsize x [1:360]

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;  % large number of colors to create the mask
    CreateMapping(E.calibrateAmp==2);

    for SS1 = 1:3
        for SS2 = 1:3
            for testedArray = 1:2
                for iii = 1:length(InterArrayInterval)

                    fdistance = zeros(E.ntrials, 1);  % distance (target, response)
                    Bstrength = zeros(E.ntrials, 2); 

                    for trial = 1:E.ntrials

                        map = struct('FX', zeros(C.nc));   % feature map
                        Map = repmat(map, C.nfeatures, 1);
                        W = CreateConnections(C.nfeatures);
                        G = zeros(1, P.nb);  % gate-closing units
                        GW = zeros(1, P.nb); % gate-closing weights

                        L = zeros(2, floor(C.nloc/2));
                        L(1,:) = randperm(floor(C.nloc/2));   % all on the first semi-circle
                        L(2,:) = randperm(floor(C.nloc/2)) + ceil(C.nloc/2);   % all on the second semi-circle % all in the second semi-circle
                        features = randperm(C.nstim);      %shuffle object features
                        F = [features(1:3); features(4:6)];

                        % first array
                        E.RI = InterArrayInterval(iii) - E.prestime;
                        [Map, W, G, GW, Focus, Afocus, content, context, Inpos, Strength, bstrength, CTime, SpatAttn] = IMencoding(Map, W, G, GW, L(1,:), F(1,:), SS1, 1);
                        usedTime = SS1*CTime; %CTime is the mean consolidation time taken
                        overTime = max(0, usedTime - InterArrayInterval(iii)); 
                        Bstrength(trial,1) = mean(bstrength);
                        
                        % second array
                        E.RI = 0.5; % check with Jacob!
                        [Map, W, G, GW, Focus, Afocus, content, context, Inpos, Strength, bstrength, CTime, SpatAttn] = IMencoding(Map, W, G, GW, L(2,:), F(2,:), SS2, 1, overTime);
                        Bstrength(trial,2) = mean(bstrength);

                        % test
                        probed = 1;  % for now
                        probestim = []; probeIdx = []; 
                        [response, rt, Map, W, G, Focus, CWcolor] = IMretrieve(Map, W, G, Focus, Afocus, probed, 1, L(testedArray,:), F(testedArray,:), probestim, probeIdx);
                        fdistance(trial) = wrap(response-F(testedArray,1), 180);   %calculate distance between response and true feature in feature space (degrees!)

                    end

                    Mdevobs(id, SS1, SS2, testedArray, iii) = mean(abs(fdistance));  %mean deviation (average over trials)
                    Mbstrength(id, SS1, SS2, :, iii) = mean(Bstrength,1)';

                    disp('    ID      Array tested   SS1      SS2     III     error');
                    disp([id, testedArray, SS1, SS2, iii, mean(Mdevobs(id, SS1, SS2, testedArray, iii))]);


                end % III
            end % testedArray
        end  % for SS2
    end % for SS1

end  % for ID

% Plot Mean(Deviation) as function of III and tested array for each combination of SS1 and SS2
PreFigure;
plotIdx = 1;
for SS1 = 1:3
    for SS2 = 1:3
        subplot(3,3,plotIdx);
        plotvector = squeeze(mean(Mdevobs(:,SS1,SS2,:,:),1));
        plot(InterArrayInterval, plotvector);
        PostFigure([-0.1, max(InterArrayInterval)+0.1, 0, 90], 'Inter-Item Interval', 'Deviation (Deg)', ['SS1=', mat2str(SS1), '; SS2=', mat2str(SS2)], {'First Array', 'Second Array'});
        plotIdx = plotIdx+1;
    end
end

% Plot Mean(Binding Strength) as function of III and tested array for each combination of SS1 and SS2
PreFigure;
plotIdx = 1;
for SS1 = 1:3
    for SS2 = 1:3
        subplot(3,3,plotIdx);
        plotvector = squeeze(mean(Mbstrength(:,SS1,SS2,:,:),1));
        plot(InterArrayInterval, plotvector);
        PostFigure([-0.1, max(InterArrayInterval)+0.1, 0, 0.2], 'Inter-Item Interval', 'B-Strength', ['SS1=', mat2str(SS1), '; SS2=', mat2str(SS2)], {'First Array', 'Second Array'});
        plotIdx = plotIdx+1;
    end
end

halt = 1;
end

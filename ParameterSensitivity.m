function [] = ParameterSensitivity(Model, ParName, Values)
% Simulation of effect of varying a parameter in a sequential-presentation
% continuous-reproduction paradigm

global P
global E
global C

setsize = 4;
E.presentation = 2;  % sequential!

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% initialize observed mean Deviation
Mdevobs = NaN(E.nsubj, length(Values), setsize, setsize);  % id, setsize, inpos, outpos
CircSD = NaN(E.nsubj, length(Values), setsize, setsize);  % id, setsize, inpos, outpos
Mrt = NaN(E.nsubj, length(Values), setsize, setsize);  % id, setsize, inpos, outpos

[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, setsize), 1:360);  %Colors = E.ntrials x E.maxsetsize x [1:360]

PP = P; % keep original Parameter structure

for parIdx = 1:length(Values)

    P = PP; 
    eval(['P.', ParName, ' = Values(parIdx);']);
    ParX = CreateIndDiff;

    for id = 1:E.nsubj

        % extract parameter values for each subject - for those parameters that vary between subjects
        for ii = 1:length(C.indVar)
            eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
        end

        % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
        CreateStimuli;
        CreateMapping(E.calibrateAmp==2);
        E.outsize = setsize;
        fdistance = NaN(E.ntrials, setsize, E.outsize);  % distance (target, response) for each input position (position in presentation order) x output position (position in test order)
        RT = zeros(1,E.ntrials, setsize, E.outsize);       % RT for each input x output position

        for trial = 1:E.ntrials

            output = Model(P, setsize, 1);  % cueing = 1 (no cue)
            for outpos = 1:E.outsize
                fdistance(trial, output.Inpos(outpos), outpos) = wrap(output.response(outpos)-output.F(outpos), 180);   %calculate distance between response and true feature in feature space (degrees!)
                RT(trial, output.Inpos(outpos), outpos) = output.rt(outpos);
            end
            tcount = tcount+1;

        end

        for inpos = 1:setsize
            for outpos = 1:E.outsize
                Mdevobs(id, parIdx, inpos, outpos) = nanmean(abs(fdistance(:,inpos,outpos)));  %mean deviation
                CircSD(id, parIdx, inpos, outpos) = circ_std(deg2rad(fdistance(~isnan(fdistance(:,inpos,outpos)),inpos,outpos)));
                Mrt(id, parIdx, inpos, outpos) = mean(RT(:,inpos,outpos))*setsize;
            end
        end

        disp('ID  Parameter  Value');
        disp([mat2str(id), '   ' ParName, '   ', mat2str(Values(parIdx))]);

    end %for ID
end  % for parIdx

% Plot Mean(Deviation) as function of set size

% Stand-alone figures of errors, averaged over oputut positions (or input
% positions for plot of output position effect)
legendtext = {'In=1', 'In=2', 'In=3','In=4', 'In=5', 'In=6', 'In=7', 'In=8'};
PreFigure;
plotvector = squeeze(nanmean(nanmean(Mdevobs,4),1));  % average over outpos (4) and subjects (1)
plot(Values, plotvector);
PostFigure([min(Values)-0.05*max(Values), 1.05*max(Values), 0, 1.05*max(max(plotvector))], ParName, 'Deviation (Deg)', [], legendtext);




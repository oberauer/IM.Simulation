function [] = SetsizeSerialRecall(Model, maxSetsize)
% Simulation of Set-size and serial recall

global P
global E
global C

E.maxsetsize = maxSetsize;
E.presentation = 2;
E.forwardrecall = 1;
%P.pgrad = P.pgrad * P.pgradBoost; % because of (expectation of) multiple tests

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

Mdevobs = NaN(E.nsubj, E.maxsetsize, E.maxsetsize);  % id, inpos, setsize
Mwact = NaN(E.nsubj, E.maxsetsize);  % id, setsize

[aa, bb, Colorgrid] = ndgrid(ones(1,E.ntrials), ones(1, E.maxsetsize), 1:360);  %Colors = E.ntrials x E.maxsetsize x [1:360]

for id = 1:E.nsubj
    
    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end
    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);
    
    for ssidx = 0:E.maxsetsize
        
        setsize = max(1, ssidx);
        E.outsize = setsize;
        fdistance = zeros(E.ntrials, setsize);  % distance (target, response) for each serial position
        wact = zeros(1,E.ntrials);
        
        for trial = 1:E.ntrials
            output = Model(P, setsize, 1);  % cueing = 1 (no cue)
            for sp = 1:setsize
                fdistance(trial, sp) = wrap(output.response(sp)-output.F(sp), 180);   %calculate distance between response and true feature in feature space (degrees!)
            end
            wact(trial) = sum(sum(abs(output.wx))); % sum of activation in weight matrix (in case of model 8, wfocus) -> CDA
        end
        for sp = 1:setsize
            Mdevobs(id, setsize, sp) = mean(abs(fdistance(:,sp)));
        end
        Mwact(id, setsize) = mean(wact);
        
    end %for setsize
    
end  % for ID

% Plot Mean(Deviation) as a function of set size and output position (simultaneous presentation)
PreFigure;
legendtext = {'SS=1', 'SS=2', 'SS=3','SS=4', 'SS=5', 'SS=6', 'SS=7', 'SS=8'};
plotvector = squeeze(mean(Mdevobs,1));
plot(plotvector');
PostFigure([0.8,setsize+0.2, 0, 1.05*max(max(plotvector))], 'Serial Position', 'Deviation (Deg)', [], legendtext);

% Plot CDA as a function of set size

CDA = mean(Mwact);
PreFigure;
ymax = max(CDA);
plot(1:E.maxsetsize, CDA);
PostFigure([0.8, E.maxsetsize+0.2, 0, 1.05*ymax], 'Set Size', 'CDA');

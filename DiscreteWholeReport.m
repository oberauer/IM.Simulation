function D = DiscreteWholeReport
% Simulation of discrete 9-AFC whole report (Adam)

global E
global C

setsize = 6;
E.test = 3;  % n-AFC
E.respAlt = [1, 3*ones(1, setsize-1), 2*ones(1, 9-setsize)]; % correct and extraset lure
E.outsize = setsize;

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

% generate parameters with individual differences
ParX = CreateIndDiff;

ProportionK = NaN(E.nsubj, setsize+1);  %

for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;  % large number of colors to create the mask
    CreateMapping(E.calibrateAmp==2);

    accuracy = zeros(E.ntrials, setsize);  % correctness

    for trial = 1:E.ntrials
        testorder = randperm(E.outpos);
        for probe = testorder
            accuracy(trial, probe) = response(1)==probe;  % probes 1:setsize are memory items F(1:setsize)
        end
    end
    NumberCorrect = sum(accuracy, 2);
    for k = 0:setsize
        ProportionK(id, k+1) = sum(NumberCorrect==k);  % proportion of trials with k correct responses
    end
    disp('    ID  0   1   2   3   4   5   6   ');
    disp([id, ProportionK]);
end  % for ID

% Plot Proportion of Number of Correct Responses
PreFigure;
plot(1:setsize, mean(ProportionK));
PostFigure([-0.5, setsize+0.5, 0, 1], 'Number Correct', 'Proportion of Trials');

D.ProportionK = ProportionK;


halt = 1;
end

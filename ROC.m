function D = ROC(Model, ptype)
% Simulation of change localization with varying number of response options
% to reconstruct ROC with the method of He, Kellen, and Singmann (2026)

global P
global E
global C

% Calibrate amplification factor on population level, if desired
if E.calibrateAmp == 1
    CreateStimuli;
    CreateMapping(1);
end

E.material = 2; 
E.test = 4; % change localization 
E.ptype = ptype; % kind of change (2 = new, 3 = swap)
setsize = 6; 
ResponseSetSize = 2:(setsize-1);

% generate parameters with individual differences
ParX = CreateIndDiff;

% Initializing some container matrices
%Choice = zeros(E.nsubj, length(ResponseSetSize));  % Choice of the changed target
PC = zeros(E.nsubj, max(ResponseSetSize));    % Proportion correct
PC(:,1) = 1; % for 1-AFC, p(correct) must be set to 1 for the ROC reconstruction
RT = zeros(E.nsubj, max(ResponseSetSize));    % Response time
CumPhit = zeros(E.nsubj, max(ResponseSetSize) + 1);
CumPFA = zeros(E.nsubj, max(ResponseSetSize) + 1);
for id = 1:E.nsubj

    % extract parameter values for each subject - for those parameters that vary between subjects
    for ii = 1:length(C.indVar)
        eval(['P.', C.indVar{ii}, ' = ParX(id, ii);']);
    end

    % for each subject, create stimuli, and an individual set of feature categories, and the corresponding mappings
    CreateStimuli;
    CreateMapping(E.calibrateAmp==2);

    for rss = 1:length(ResponseSetSize)

        E.rss = ResponseSetSize(rss);
        choice = zeros(E.ntrials, 1);
        rt = zeros(E.ntrials, 1);

        for trial = 1:E.ntrials
            output = Model(P, setsize, 1);   % run model on 1 trial, returns predictions (output is a structure with lots of variables in it)
            choice(trial) = output.response;  % the first entry of response is the actual response
            rt(trial) = output.rt;    % response time
        end

        PC(id, ResponseSetSize(rss)) = mean(choice==1);
        RT(id, ResponseSetSize(rss)) = mean(rt);

        disp(['      ID      RSS       PC        RT']);
        disp([id, E.rss, PC(id, ResponseSetSize(rss)), RT(id, ResponseSetSize(rss))]);

    end

    %%% compute probability that the changed item is ranked i (Eq. 2)
    % fit Block-Marshak Inequalities to PC (advice by Y. He and David Kellen, May 28, 2026)
    PC_BM = ones(1,5);
    PC_BM(2:5) = BlockMarshakFit(PC(id,2:5));
    m = max(ResponseSetSize);
    R = zeros(1,m);
    for i = 1:m
        Sum = 0; 
        for j = (m-i+1):m
            Sum = Sum + nchoosek(i-1, j-(m-i+1)) * (-1).^(j-(m-i+1)) * PC_BM(j);
        end
        R(i) = nchoosek(m-1, i-1) * Sum;
    end
    % compute probability that a lure is ranked i (Eq. 4)
    Q = (1-R)./(m-1);
    % construct ROC
    CumPhit(id,:) = cumsum([0, R]); 
    CumPFA(id,:) = cumsum([0, Q]);

end  % for ID

grayColor = [.7 .7 .7];
PreFigure;
hold on
for id = 1:E.nsubj
    plot(CumPFA(id,:), CumPhit(id,:), 'color', grayColor); 
end
plot(mean(CumPFA), mean(CumPhit), 'color', 'red');
PostFigure([0, 1, 0, 1], 'P(FA)', 'P(Hit)');

PreFigure;
hold on
for id = 1:E.nsubj
    plot(PC(id,:), 'color', grayColor); 
end
plot([1, ResponseSetSize], mean(PC), 'color', 'red');
PostFigure([0, max(ResponseSetSize)+1, 0, 1.1], 'RSS', 'P(Correct)');

D.CumPhit = CumPhit;
D.CumPFA = CumPFA;
D.PC = PC;
D.RT = RT; 
output = 1; 


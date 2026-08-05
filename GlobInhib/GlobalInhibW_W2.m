%%% Self-activation and global inhibition in W;
% Explore effects of encTime and duration on SPC
% With a second W2 that gradually learns the first

% Reducing inhibition --> more primacy, better performance
% Increasing enctime = duration: less primacy, more recency
% Increasing input drive: At short encTime it slightly improves recency w/o
% hurting primacy. At long encTime it hurts primacy with little effect on
% recency
% Increasing selfAct does little (slightly decreases performance)
% Increasing inhibition reduces primacy, improves recency, especially at
% longer encTimes
% Higher asymptote tilts SPCs towards more recency, especially at longer
% encTimes
% Lower encRate decreases recency

clear variables
%close all;

nTrials = 5000;
setsize = 6;
EncTime = 0.5;
Duration = [0.1, 0.2, 0.3, 0.4, 0.8];
tstep = 0.02;
shunting = 1;

inputDrive = 1.5;
kappa = 25;
kappaCat = 10;
strengthSD = 0.1;
dnoise = 0.01;
nCat = 8;
encRate = 100;
lrate = 0.0;
primacy = 0.2; 
commit = 0;   % strength of commitment --> how much inhibition of the committed units is dampened
cThreshold = 0.01;

if shunting == 0
    selfAct = 1;
    inhib = 1;
    asyW = 3;
end
if shunting == 0.5
    selfAct = 1;
    inhib = 0.5;
    asyW = 3;
end
if shunting == 1
    selfAct = 1.0;   % hand-set
    inhib = 0.7;   % hand-set
    asyW = 4;
end

stepSize = round(360/nCat);
xradCat = deg2rad(stepSize:stepSize:360);
xrad = deg2rad(1:360);

WCat = zeros(360, nCat);
for cat = 1:nCat
    WCat(:,cat) = VonMisesN(xrad, deg2rad(cat*stepSize), kappaCat)';
end

SPC = NaN(length(Duration), setsize);
meanError = NaN(1, length(Duration));

for dIdx = 1:length(Duration)
    duration = Duration(dIdx);
    encTime = min(duration, EncTime);
    nSteps = round(duration./tstep);

    Error = zeros(nTrials, setsize);
    maxW = zeros(nTrials, nSteps);

    for trial = 1:nTrials
        W = zeros(nCat);
        W2 = zeros(nCat);
        Committed = zeros(nCat);
        InputW = zeros(nCat);
        stim = randperm(360, setsize);
        loc = randperm(360, setsize);
        strength = max(0, exp(-primacy*(0:(setsize-1))) + randn(1, setsize) * strengthSD);

        %  Encoding
        for item = 1:setsize

            Stim = VonMises(xrad, deg2rad(stim(item)), kappa) * WCat;
            Loc = VonMises(xrad, deg2rad(loc(item)), kappa) * WCat;
            InputW = strength(item) * Loc' * Stim;

            if shunting == 0.0, GI = @(x,t,Input,C) selfAct.*x + (t <= encTime/tstep)*inputDrive*(asyW-max(x(:)))*Input - (1-C)*inhib*sum(x(:)); end % non-shunting version
            if shunting == 0.5, GI = @(x,t,Input,C) selfAct.*x.*(asyW-x) + (t <= encTime/tstep)*inputDrive*(asyW-max(x(:)))*Input - (1-C)*inhib*sum(x(:)); end % half-shunting version
            if shunting == 1.0, GI = @(x,t,Input,C) selfAct.*x.*(asyW-x) + (t <= encTime/tstep)*inputDrive*(asyW-max(x(:)))*Input - x.*(1-C)*inhib*sum(x(:)); end  %shunting version

            for t = 1:nSteps
                InW = (1 - exp(-encRate*t*tstep)) * InputW;

                % RK4
                k1 = GI(W,t,InputW,Committed);
                k2 = GI(max(0, min(asyW, W + 0.5*tstep*k1)),t,InputW,Committed);
                k3 = GI(max(0, min(asyW, W + 0.5*tstep*k2)),t,InputW,Committed);
                k4 = GI(max(0, min(asyW, W + tstep*k3)),t,InputW,Committed);
                W = max(0, min(asyW, W + (tstep/6) * (k1 + 2*k2 + 2*k3 + k4)));
                maxW(trial, t) = max(W(:));
                W2 = W2 + lrate*W;
            end

            Committed(W > cThreshold) = commit;

        end

        % Retrieval
        for item = 1:setsize

            cue = VonMises(xrad, deg2rad(loc(item)), kappa) * WCat;
            reAct = cue * (W+W2) * WCat';
            reActN = reAct + randn(1,360)*dnoise;
            retrieved = find(reActN == max(reActN));
            Error(trial, item) = abs(wrap(retrieved - stim(item), 180));

        end

    end

    SPC(dIdx, 1:setsize) = mean(Error);
    meanError(dIdx) = mean(Error(:));

end

PreFigure([], [], 2);
subplot(1,2,1);
plot(1:setsize, SPC');
PostFigure([0.5, setsize+0.5, 0, 90], 'Input Position', 'Error', [], vec2legend(Duration));
subplot(1,2,2);
plot(Duration, meanError);
PostFigure([0, max(Duration)+0.1, 0, 90], 'Duration', 'Error');

% MaxW = mean(maxW);
% Time = tstep*(1:nSteps);
% subplot(2,3,setsize);
% plot(Time, MaxW);
% PostFigure([0, duration, 0, max(0.1, max(MaxW))], 'Time (s)', 'max(W)', num2str(setsize) );











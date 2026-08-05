%%% Self-activation and global inhibition in W;
% with a binding layer B


clear variables
%close all;

nTrials = 5000;
setsize = 6;
EncTime = 0.1;
Duration = [0.1, 0.2, 0.3, 0.4, 0.8];
tstep = 0.02;
shunting = 1;

inputDrive = 3;
kappa = 25;
kappaB = 25; 
kappaCat = 10;
strengthSD = 0.1;
dnoise = 0.01;
nCat = 8;
nb = 90; 
encRate = 10;
primacy = 0.01;  % rate of exponential decline of encoding strength (smaller -> flatter)
commit = 0.75;   % strength of commitment --> how much inhibition of the committed units is dampened
cThreshold = 0.1;  %threshold of activation in W: Only units above the threshold are committed

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
    selfAct = 1.0;   
    inhib = 1;  
    asyW = 4;
end

stepSize = round(360/nCat);
xradCat = deg2rad(stepSize:stepSize:360);
xrad = deg2rad(1:360);
stepSizeB = round(360/nb);
xradB = deg2rad(stepSizeB:stepSizeB:360);
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
        W = zeros(2*nCat, nb);
        InputW = zeros(2*nCat, nb);
        Committed = zeros(2*nCat, nb);
        stim = randperm(360, setsize);
        loc = randperm(360, setsize);
        strength = max(0, exp(-primacy*(0:(setsize-1))) + randn(1, setsize) * strengthSD);

        %  Encoding
        for item = 1:setsize

            Stim = VonMises(xrad, deg2rad(stim(item)), kappa) * WCat;
            Loc = VonMises(xrad, deg2rad(loc(item)), kappa) * WCat;
            B = VonMises(xradB, randperm(360,1), kappaB);
            InputW = strength(item) * [Loc, Stim]' * B;

            if shunting == 0.0, GI = @(x,t,Input,C) selfAct.*x + (t <= encTime/tstep)*inputDrive*(asyW-max(x(:)))*Input - (1-C)*inhib*sum(x(:)); end % non-shunting version
            if shunting == 0.5, GI = @(x,t,Input,C) selfAct.*x.*(asyW-x) + (t <= encTime/tstep)*inputDrive*(asyW-max(x(:)))*Input - (1-C)*inhib*sum(x(:)); end % half-shunting version
            if shunting == 1.0, GI = @(x,t,Input,C) selfAct.*x.*(asyW-x) + (t <= encTime/tstep)*inputDrive*(asyW-max(x(:)))*Input - x.*(1-C)*inhib*sum(x(:)); end  %shunting version

            for t = 1:nSteps
                InW = (1 - exp(-encRate*t*tstep)) * InputW;

                % RK4
                k1 = GI(W,t,InW,Committed);
                k2 = GI(max(0, min(asyW, W + 0.5*tstep*k1)),t,InW,Committed);
                k3 = GI(max(0, min(asyW, W + 0.5*tstep*k2)),t,InW,Committed);
                k4 = GI(max(0, min(asyW, W + tstep*k3)),t,InW,Committed);
                W = max(0, min(asyW, W + (tstep/6) * (k1 + 2*k2 + 2*k3 + k4)));
                maxW(trial, t) = max(W(:));
            end

            Committed(W > cThreshold) = commit;

        end

        % Retrieval
        for item = 1:setsize

            cueLoc = VonMises(xrad, deg2rad(loc(item)), kappa) * WCat;
            cue = [cueLoc, zeros(1, nCat)];
            retrievedB = cue*W;
            retrievedVec = retrievedB*W';
            reAct = retrievedVec((nCat+1):end) * WCat';
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











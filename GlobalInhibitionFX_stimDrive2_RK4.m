%%% Try out self-activation and global inhibition in FX to implement decay
% Variant including the stimulus drive in the differential equation
% This is using Runge-Kutta 4

clear variables
%close all;

nTrials = 30;
maxSetsize = 12;
stimDrive = 5;
encTime = 0.15;
tstep = 0.02;
duration = 1; % in seconds
nSteps = round(duration./tstep);
shunting = 0;

if shunting == 0
    selfAct = 1.0;
    inhib = 0.002;
    asyFX = 5;
end
if shunting == 0.5
    selfAct = 1;
    inhib = 0.002;
    asyFX = 5;
end
if shunting == 1
    selfAct = 5;   % hand-set
    inhib = 0.0003;   % hand-set
    asyFX = 5;
end
kappa = 25;
strengthSD = 0.2;
threshold = 0.1;
xrad = deg2rad(1:360);
NumAlive = zeros(nSteps, maxSetsize);
MaxAct = zeros(nSteps, maxSetsize);
MeanAct = zeros(nSteps, maxSetsize);
SumAct = zeros(nSteps, maxSetsize);

tic
for setsize = 1:maxSetsize

    numberAlive = zeros(nTrials, nSteps);
    maxAct = zeros(nTrials, nSteps);
    meanPeakAct = zeros(nTrials, nSteps);
    sumAct = zeros(nTrials, nSteps);

    for trial = 1:nTrials
        FX = zeros(360);
        Input = zeros(360);
        stim = randperm(360, setsize);
        loc = randperm(360, setsize);
        strength = max(0, randn(1, setsize) * strengthSD + 1);

        % stimulus input
        for item = 1:setsize
            Input = Input + strength(item) * VonMises(xrad, deg2rad(loc(item)), kappa)' * VonMises(xrad, deg2rad(stim(item)), kappa);
        end


        % dynamics
        % if trial == 1 && setsize == maxSetsize
        %     PreFigure;
        %     plotIdx = 1;
        % end

        if shunting == 0.0, GI = @(x,t) selfAct.*x + (t <= encTime/tstep)*stimDrive*(asyFX-max(x(:)))*Input - inhib*sum(x(:)); end % non-shunting version
        if shunting == 0.5, GI = @(x,t) selfAct.*x.*(asyFX-x) - inhib*sum(x(:)); end % half-shunting version
        if shunting == 1.0, GI = @(x,t) selfAct.*x.*(asyFX-x) - x.*inhib*sum(x(:)); end  %shunting version

        for t = 1:nSteps

            maxAct(trial, t) = max(FX(:));
            alive = 0;
            sumPeakact = 0;
            for item = 1:setsize
                sumPeakact = sumPeakact + FX(loc(item), stim(item));
                alive = alive + round(FX(loc(item), stim(item)) > threshold);
            end
            meanPeakAct(trial, t) = sumPeakact./setsize;
            sumAct(trial,t) = sum(FX(:));

            % RK4
            k1 = GI(FX,t);
            k2 = GI(max(0, min(asyFX, FX + 0.5*tstep*k1)),t);
            k3 = GI(max(0, min(asyFX, FX + 0.5*tstep*k2)),t);
            k4 = GI(max(0, min(asyFX, FX + tstep*k3)),t);
            FX = max(0, min(asyFX, FX + (tstep/6) * (k1 + 2*k2 + 2*k3 + k4)));

            numberAlive(trial, t) = alive;
            % if setsize == maxSetsize
            %     if trial == 1 && mod(t, 10) == 0 && plotIdx < 13
            %         subplot(3,4,plotIdx);
            %         image(FX*100);
            %         plotIdx = plotIdx + 1;
            %     end
            % end
        end
    end
    NumAlive(:, setsize) = mean(numberAlive, 1)';
    MaxAct(:, setsize) = mean(maxAct, 1)';
    MeanAct(:, setsize) = mean(meanPeakAct, 1)';
    SumAct(:, setsize) = mean(sumAct, 1)';
end
toc

Time = tstep*(1:nSteps);
PreFigure([], [], 2);
subplot(2,2,1);
plot(Time, NumAlive);
PostFigure([0, max(Time), 0, setsize+1], 'Time (s)', 'Number of Items Alive', [], vec2legend(1:maxSetsize));
subplot(2,2,3);
plot(Time, MaxAct);
PostFigure([0, max(Time), 0, max(0.1, max(MaxAct(:)))], 'Time (s)', 'Max. Peak Act.');
subplot(2,2,4);
plot(Time, MeanAct);
PostFigure([0, max(Time), 0, max(0.1, max(MeanAct(:)))], 'Time (s)', 'Mean Peak Act.');

subplot(2,2,2);
plot(Time, SumAct);
PostFigure([0, max(Time), 0, max(0.1, max(SumAct(:)))], 'Time (s)', 'Sum Act.');

% SumAct1 = SumAct(round(nSteps/2), :);
% subplot(2,2,2);
% plot(1:maxSetsize, SumAct1);
% PostFigure([0, maxSetsize+1, 0, 1.2*max(SumAct1(:))], 'Set Size', 'Summed Act.');





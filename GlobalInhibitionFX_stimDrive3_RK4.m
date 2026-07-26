%%% Try out self-activation and global inhibition in FX to implement decay
% Variant including the stimulus drive in the differential equation
% Variant also including IRO in the differential equation
% This is using Runge-Kutta 4

clear variables
close all;

nTrials = 30;
maxSetsize = 12;
stimDrive = 3;
encTime = 0.15;
tstep = 0.02;
duration = 1; % in seconds
ior = 0.5;    % inhibition of return

%ior = 0; 

nSteps = round(duration./tstep);
shunting = 0;

if shunting == 0
    selfAct = 2;
    inhib = 0.004;
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
MeanPeakAct = zeros(nSteps, maxSetsize);
MeanAct = zeros(nSteps, maxSetsize);
Alpha = zeros(nSteps, maxSetsize);
eW = 0.3 + randn(360, 52);
eNoise = 0.25;


% 
% selfAct = selfAct*4;
% inhib = inhib * 4;
% stimDrive = stimDrive*4;

tic
for setsize = 1:maxSetsize

    numberAlive = zeros(nTrials, nSteps);
    maxAct = zeros(nTrials, nSteps);
    meanPeakAct = zeros(nTrials, nSteps);
    sAttn = zeros(nTrials, nSteps);
    alpha = zeros(nTrials, nSteps);
    peak1 = zeros(nTrials, nSteps);
    peak2 = zeros(nTrials, nSteps);

    for trial = 1:nTrials
        FX = zeros(360);
        Input = zeros(360);
        IORinput = zeros(360);
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

        if shunting == 0.0, GI = @(x,t,ior) selfAct.*x + (t <= encTime/tstep)*stimDrive*(asyFX-max(x(:)))*Input - inhib*sum(x(:)) - stimDrive*ior; end % non-shunting version
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
            spatAttn = mean(FX,2);
            sAttn(trial,t) = sum(spatAttn);
            alpha(trial,t) = sum(abs(spatAttn' * eW + randn(1,52)*eNoise));
            if setsize == 2
                peak1(trial,t) = FX(loc(1), stim(1));
                peak2(trial,t) = FX(loc(2), stim(2)); 
            end

            % inhibition of return: one item every 300 ms, in order of
            % presentation, lasting for 300 ms (until IOR of the next item
            % is initiated)
            for item = 1:setsize
                if t == round((item*0.3)/tstep)
                    IORinput = ior * (VonMises(xrad, deg2rad(loc(item)), kappa)' * VonMises(xrad, deg2rad(stim(item)), kappa)); 
                end
            end

            % RK4
            k1 = GI(FX,t,IORinput);
            k2 = GI(max(0, min(asyFX, FX + 0.5*tstep*k1)),t,IORinput);
            k3 = GI(max(0, min(asyFX, FX + 0.5*tstep*k2)),t,IORinput);
            k4 = GI(max(0, min(asyFX, FX + tstep*k3)),t,IORinput);
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
    MeanPeakAct(:, setsize) = mean(meanPeakAct, 1)';
    MeanAct(:, setsize) = mean(sAttn, 1)';
    Alpha(:, setsize) = mean(alpha, 1)';
    if setsize == 2
        Peaks(:, 1) = mean(peak1);
        Peaks(:, 2) = mean(peak2);
    end
    
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
plot(Time, MeanPeakAct);
PostFigure([0, max(Time), 0, max(0.1, max(MeanPeakAct(:)))], 'Time (s)', 'Mean Peak Act.');

subplot(2,2,2);
plot(Time, MeanAct);
PostFigure([0, max(Time), 0, max(0.1, max(MeanAct(:)))], 'Time (s)', 'Mean Act.');

% SumAct1 = SumAct(round(nSteps/2), :);
% subplot(2,2,2);
% plot(1:maxSetsize, SumAct1);
% PostFigure([0, maxSetsize+1, 0, 1.2*max(SumAct1(:))], 'Set Size', 'Summed Act.');


Alpha400 = mean(Alpha(round(0.4/tstep):round(1/tstep),:));  % 400 to 1000 ms
PreFigure([], [], 2);
subplot(1,2,1);
plot(Time, Alpha);
PostFigure([0, max(Time), 0, max(0.1, max(Alpha(:)))], 'Time (s)', 'Alpha');
subplot(1,2,2);
plot(1:setsize, Alpha400);
PostFigure([0, maxSetsize+1, 0, max(0.1, 0.3*max(Alpha(:)))], 'Set Size', 'Alpha(400-1000)');

PreFigure([], [], 2);
plot(Time, Peaks');
PostFigure([0, max(Time), 0, 1.1*max(Peaks(:))], 'Time (s)', 'Peak Act.', 'Set Size 2', {'Item 1', 'Item 2'});
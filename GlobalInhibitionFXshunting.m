%%% Try out self-activation and global inhibition in FX to implement decay
% Variant with shunting equations: inhibition is proportional to currenct
% value 

clear variables
close all;

nTrials = 100;
maxSetsize = 12;
%tstep = 0.05;  %50 ms steps
tstep = 0.01;

if tstep == 0.01
    nSteps = 200; % 10-millisecond steps
    selfAct1 = 1 + 1.5*tstep; % general rule for any tstep
    inhib1 = 0.004*tstep;     % general rule
    selfAct = 0.015;   % hand-set
    inhib = 0.00004;   % hand-set
end
if tstep == 0.05
    nSteps = 40; % 50-millisecond steps
    selfAct1 = 1 + 1.5*tstep; % general rule for any tstep (but does not work)
    inhib1 = 0.004*tstep;     % general rule
    selfActClaude = 1.077;          % recommendation by Claude
    inhibClaude = 0.0002;           % recommendation by Claude
    selfAct = 0.03;
    inhib = 0.00008;
end

kappa = 25;
strengthSD = 0.1;
asyFX = 4; 
threshold = 0.1;
xrad = deg2rad(1:360);

NumAlive = zeros(nSteps, maxSetsize);
MaxAct = zeros(nSteps, maxSetsize);
MeanAct = zeros(nSteps, maxSetsize);

for setsize = 1:maxSetsize

    numberAlive = zeros(nTrials, nSteps);
    maxAct = zeros(nTrials, nSteps);
    meanAct = zeros(nTrials, nSteps);

    for trial = 1:nTrials
        FX = zeros(360);
        stim = randperm(360, setsize);
        loc = randperm(360, setsize);
        strength = max(0, randn(1, setsize) * strengthSD + 1);
        % encoding
        for item = 1:setsize
            FX = FX + strength(item) * VonMises(xrad, deg2rad(loc(item)), kappa)' * VonMises(xrad, deg2rad(stim(item)), kappa);
        end
        % dynamics
        % if trial == 1 && setsize == maxSetsize
        %     PreFigure;
        %     plotIdx = 1;
        % end
        for t = 1:nSteps
            maxAct(trial, t) = max(FX(:));
            alive = 0;
            meanact = 0;
            for item = 1:setsize
                meanact = meanact + FX(loc(item), stim(item));
                alive = alive + round(FX(loc(item), stim(item)) > threshold);
            end
            meanAct(trial, t) = meanact./setsize;
            FX = min(asyFX, FX + (asyFX-FX).*selfAct.*FX);
            summedAct = sum(FX(:));
            FX = max(0, FX - FX.*(inhib*summedAct)); 
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
    MeanAct(:, setsize) = mean(meanAct, 1)';
end

Time = tstep*(1:nSteps);
PreFigure([], [], 2);
subplot(2,2,1);
plot(Time, NumAlive);
PostFigure([0, max(Time), 0, setsize+1], 'Time (s)', 'Number of Items Alive', [], vec2legend(1:maxSetsize));
subplot(2,2,3);
plot(Time, MaxAct);
PostFigure([0, max(Time), 0, max(MaxAct(:))], 'Time (s)', 'Max. Activation',  [], vec2legend(1:maxSetsize));
subplot(2,2,4);
plot(Time, MeanAct);
PostFigure([0, max(Time), 0, max(MeanAct(:))], 'Time (s)', 'Mean Activation',  [], vec2legend(1:maxSetsize));
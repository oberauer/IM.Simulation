%%% How do dnoise and distractors affect behavior of SDM?

clear all
close all

kappa = 10;
strengthening = 0.5;
Boundary = 30:10:70;
Dnoise = 0.2:0.1:0.8;
scaling = 1; % scaling factor to simulate the reduction of binding resource

nSubj = 30;
nTrials = 200;
setsize = 2;
x = pi*(-179:180)./180;
distshift = 100;   % shift (in degrees) of the distractor relative to the target (the target is always at 180)
diststrength = 0.75;
ntstrength = 0.75;
tstep = 0.05;
nsteps = 5000*tstep;   % that should be more than enough (5 sec)
option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

for dIdx = 1:length(Dnoise)
    for dist = 1:2
        for cue = 1:2
            dnoise = Dnoise(dIdx);
            ntfeature = zeros(nTrials, setsize);
            A = zeros(nTrials, 360);
            P = zeros(nTrials, 360);
            Pt = zeros(nTrials, 1); 
            Pd = zeros(nTrials, 1); 
            Pnt = zeros(nTrials, 1); 
            for trial = 1:nTrials
                A(trial, :) = (1 + (cue-1)*strengthening) .* VonMises(x, 0, kappa);
                ntfeature(trial, :) = randperm(360, setsize);
                for s = 2:setsize
                    A(trial, :) = A(trial, :) + ntstrength .* VonMises(x, x(ntfeature(trial, s-1)), kappa);
                end
                A(trial, :) = A(trial, :) + (dist-1)*diststrength*VonMises(x, x(180+distshift), kappa);    % distractor after the retro-cue
                P(trial, :) = exp(A(trial,:)./dnoise)./sum(exp(A(trial,:)./dnoise));
                Pt(trial) = sum(P(trial, 170:190));
                Pd(trial) = sum(P(trial, (170:190)+distshift));
                for (s = 2:setsize), Pnt(trial) = Pnt(trial) + sum(P(trial, 181+wrap((ntfeature(trial, s-1)-10-181):(ntfeature(trial, s-1)+10-181), 180))); end
            end
            PT(dIdx, dist, cue) = mean(Pt);
            PD(dIdx, dist, cue) = mean(Pd);
            PNT(dIdx, dist, cue) = mean(Pnt);
        end
    end
end

PreFigure;
subplot(3,2,1);
plot(Dnoise, squeeze(PT(:,:,1)));
PostFigure([0, max(Dnoise), 0, 1], 'Dnoise', 'PT', 'No Cue', {'No Distr', 'Distr'}); 
subplot(3,2,2);
plot(Dnoise, squeeze(PT(:,:,2)));
PostFigure([0, max(Dnoise), 0, 1], 'Dnoise', 'PT', 'Cue', {'No Distr', 'Distr'}); 

subplot(3,2,3);
plot(Dnoise, squeeze(PD(:,:,1)));
PostFigure([0, max(Dnoise), 0, 1], 'Dnoise', 'PD', 'No Cue', {'No Distr', 'Distr'}); 
subplot(3,2,4);
plot(Dnoise, squeeze(PD(:,:,2)));
PostFigure([0, max(Dnoise), 0, 1], 'Dnoise', 'PD', 'Cue', {'No Distr', 'Distr'}); 

subplot(3,2,5);
plot(Dnoise, squeeze(PNT(:,:,1)));
PostFigure([0, max(Dnoise), 0, 1], 'Dnoise', 'PNT', 'No Cue', {'No Distr', 'Distr'}); 
subplot(3,2,6);
plot(Dnoise, squeeze(PNT(:,:,2)));
PostFigure([0, max(Dnoise), 0, 1], 'Dnoise', 'PNT', 'Cue', {'No Distr', 'Distr'}); 

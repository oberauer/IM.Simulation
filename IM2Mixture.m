%%% How do parameters of the IM simulation translate into parameters of the
%%% 3-parameter mixture model, with and without colorwheel interference and headstart for retrieval?

clear all
close all

kappa = 7;
strengthening = 0; % 0.25;  % strengthening by the cue
Boundary = 10:20:70;
Dnoise = 1:1:3;
meanScaling = 1; % scaling factor to simulate the reduction of binding resource
sdScaling = 0.8; 

nSubj = 30;
nTrials = 200;
setsize = 5;
x = pi*(-179:180)./180;
distshift = 20 * (randi(10,nTrials,1)-5);   % shift (in degrees) of the distractor relative to the target (the target is always at 180)
diststrength = 0.3;
ntstrength = 0.35;
tstep = 0.05;
nsteps = 5000*tstep;   % that should be more than enough (5 sec)
option = optimset('Display','off','TolFun',1e-10, 'FunValCheck','on', 'MaxIter', 2000);

Data.D = repmat(0:(setsize-1), nTrials, 1);
Data.Dcang = zeros(nTrials, setsize, 360);
ntshift = zeros(nTrials, setsize-1);
NTx = zeros(nTrials, setsize-1);
for j = 1:nTrials
    Data.Dcang(j, 1, :) = x;
    ntshift(j, :) = randperm(360, setsize-1);
    for s = 2:setsize
        Data.Dcang(j, s, :) = circshift(x, ntshift(j, s-1));
        NTx(j, s-1) = x(Data.Dcang(j, s, :)==0);
    end
    Data.Dwang(j,:) = circshift(x, distshift(j));
end
Data.setsize = setsize*ones(nTrials,1);

RT = zeros(nSubj, length(Boundary), length(Dnoise), 2, 2);
MMSD = zeros(nSubj, length(Boundary), length(Dnoise), 2, 2);
MMguessing = zeros(nSubj, length(Boundary), length(Dnoise), 2, 2);
MMtranspos = zeros(nSubj, length(Boundary), length(Dnoise), 2, 2);
MMcwattraction = zeros(nSubj, length(Boundary), length(Dnoise), 2, 2);


for boundIdx = 1:length(Boundary)
    for dnoiseIdx = 1:length(Dnoise)
        for distr = 1:2   % addition of a distractor (after the retro-cue, if any -> color-wheel distractor)
            for cue = 1:2   % no retro-cue vs. retro-cue (head-start of retrieval)

                response = zeros(1, nTrials);
                rt = zeros(1, nTrials);
                tic
                for id = 1:nSubj
                    scaling = max(0.1, meanScaling + randn(1,nTrials)*sdScaling); 
                    for trial = 1:nTrials
                        Adrift = (1 + (cue-1)*strengthening) .* VonMises(x, 0, kappa);
                        for s = 2:setsize, Adrift = Adrift + ntstrength .* VonMises(x, NTx(trial, s-1), kappa); end
                        Adrift = scaling(trial)*Adrift;
                        A = zeros(1, 360);
                        if (cue == 2)
                            A = sum(Adrift + randn(1000*tstep, 360)*Dnoise(dnoiseIdx));  % head start
                        end
                        if distr == 2
                            Adrift = Adrift + meanScaling*diststrength*VonMises(x, x(180+distshift(trial)), kappa);    % distractor after the retro-cue
                        end
                        A = A + cumsum(Adrift + randn(nsteps, 360)*Dnoise(dnoiseIdx));  % continued sampling
                        maxA = max(A, [], 2); % maximum value in each row of A = maximum after each time step
                        t = find(maxA > Boundary(boundIdx), 1);
                        if isempty(t), t=nsteps; end
                        response(trial) = find(A(t,:)==max(A(t,:)), 1);  % responses are ordered by output position ("probed" is incremented from 1 to E.outsize)
                        rt(trial) = t;
                    end
                    Data.response = response';
                    startparms = [15, .1, .1, .1];
                    lb = [eps, 0, 0, 0]; ub = [90, 1, 1, 1];
                    npar = 4;         %2 -> Zhang-Luck mixture, 3 -> Bays mixture, 4 -> Souza & Oberauer mixture (including colorwheel attraction)
                    MMloglik = 500000;
                    itercount = 0;
                    while MMloglik > 400000
                        [MMparms, MMloglik] = fminsearchbnd(@(x) MM(x, Data, 2), startparms(1:npar), lb(1:npar), ub(1:npar), option);
                        itercount = itercount + 1;
                    end
                    RT(id, boundIdx, dnoiseIdx, distr, cue) = mean(rt)*tstep;
                    MMSD(id, boundIdx, dnoiseIdx, distr, cue) = MMparms(1);
                    MMguessing(id, boundIdx, dnoiseIdx, distr, cue) = MMparms(2);
                    MMtranspos(id, boundIdx, dnoiseIdx, distr, cue) = MMparms(3);
                    MMcwattraction(id, boundIdx, dnoiseIdx, distr, cue) = MMparms(4);
                end
                toc
                disp('    Boundary   Dnoise   distractor cueing    LL/1000   SD       Pmem      Ptrans    Pwheel ');
                disp([boundIdx, dnoiseIdx, distr, cue, MMloglik/1000, MMparms(1), 1-sum(MMparms(2:4)), MMparms(3), MMparms(4)]);

            end
        end
    end
end

MMSD = squeeze(mean(MMSD, 1));
MMguessing = squeeze(mean(MMguessing, 1));
MMtranspos = squeeze(mean(MMtranspos, 1));
MMcwattraction = squeeze(mean(MMcwattraction, 1));

MMPm = 1 - MMtranspos - MMguessing - MMcwattraction;
DistrEffectPm = squeeze(MMPm(:,:,1,:) - MMPm(:,:,2,:)); 

PreFigure;
plotIdx = 1;
for distr = 1:2
    for cue = 1:2
        subplot(2,2,plotIdx);
        plot(Boundary, MMSD(:,:,distr,cue));
        PostFigure([0, max(Boundary)+5, 0, max(MMSD(:))+5], 'Boundary', 'Mean SD', ['Distr:', mat2str(distr), '; Cueing:', mat2str(cue)], vec2legend(Dnoise));
        plotIdx = plotIdx + 1; 
    end
end

PreFigure;
plotIdx = 1;
for distr = 1:2
    for cue = 1:2
        subplot(2,2,plotIdx);
        plot(Boundary, MMPm(:,:,distr,cue));
        PostFigure([0, max(Boundary)+5, 0, 1], 'Boundary', 'P(mem)', ['Distr:', mat2str(distr), '; Cueing:', mat2str(cue)], vec2legend(Dnoise));
                plotIdx = plotIdx + 1; 
    end
end

PreFigure;
plotIdx = 1;
for distr = 1:2
    for cue = 1:2
        subplot(2,2,plotIdx);
        plot(Boundary, MMtranspos(:,:,distr,cue));
        PostFigure([0, max(Boundary)+5, 0, 1], 'Boundary', 'P(transpos)', ['Distr:', mat2str(distr), '; Cueing:', mat2str(cue)], vec2legend(Dnoise));
                plotIdx = plotIdx + 1; 
    end
end

PreFigure;
plotIdx = 1;
for distr = 1:2
    for cue = 1:2
        subplot(2,2,plotIdx);
        plot(Boundary, MMguessing(:,:,distr,cue));
        PostFigure([0, max(Boundary)+5, 0, 1], 'Boundary', 'P(guess)', ['Distr:', mat2str(distr), '; Cueing:', mat2str(cue)], vec2legend(Dnoise));
                plotIdx = plotIdx + 1; 
    end
end

PreFigure;
plotIdx = 1;
for distr = 1:2
    for cue = 1:2
        subplot(2,2,plotIdx);
        plot(Boundary, MMcwattraction(:,:,distr,cue));
        PostFigure( [0, max(Boundary)+5, 0, 1], 'Boundary', 'P(CW)', ['Distr:', mat2str(distr), '; Cueing:', mat2str(cue)], vec2legend(Dnoise));
                plotIdx = plotIdx + 1; 
    end
end

PreFigure;
plotIdx = 1;
for distr = 1:2
    for cue = 1:2
        subplot(2,2,plotIdx);
        plot(Boundary, RT(:,:,distr,cue));
        PostFigure([0, max(Boundary)+5, 0, max(RT(:))+0.5], 'Boundary', 'Mean RT', ['Distr:', mat2str(distr), '; Cueing:', mat2str(cue)], vec2legend(Dnoise));
        plotIdx = plotIdx + 1; 
    end
end

PreFigure;
subplot(1,2,1); 
plot(Boundary, DistrEffectPm(:,:,1));
PostFigure( [0, max(Boundary)+5, -0.2, 0.4], 'Boundary', 'Distractor Effect on Pmem', 'No Cueing', vec2legend(Dnoise));
subplot(1,2,2); 
plot(Boundary, DistrEffectPm(:,:,2));
PostFigure( [0, max(Boundary)+5, -0.2, 0.4], 'Boundary', 'Distractor Effect on Pmem', 'Cueing', vec2legend(Dnoise));

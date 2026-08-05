%%% Self-activation and global inhibition in W;
% with a binding layer B
%
% ===================================================================
% CHANGES FROM ORIGINAL (see chat for full derivation/evidence):
%
% 1) commit: 0.85 -> 0.70
%    `commit` is the parameter that sets the primacy/recency balance.
%    It dampens a committed unit's OWN exposure to inhibition (via the
%    (1-C) factor) but does NOT reduce that unit's contribution to the
%    shared sum(x(:)) term that inhibits every other, not-yet-committed
%    unit. So whichever item locks in first taxes everything encoded
%    after it, and how strongly depends on `commit`:
%       commit ~0-0.3  -> recency only (last item wins)
%       commit ~0.85+  -> primacy only, later items can be driven to
%                         near-zero commitment (worst at duration=0.8,
%                         see point 2)
%       commit ~0.65-0.7 -> genuine primacy+recency U-shape, but see
%                         the caveat below: reliable mainly for
%                         Duration in ~0.2-0.4s, not the full range.
%
% 2) Added tracking of NewCommit(dIdx, item): mean number of NEWLY
%    committed W-units at each item, so you can see directly whether
%    duration is producing sparser commitment for early items, and
%    whether that's costing performance (compare against SPC).
%
% 3) NOT changed, but flagged as important: EncTime = 0.5 sits inside
%    your Duration range [0.1 0.2 0.3 0.4 0.8]. The "commit fewer
%    weights with more duration" behavior in this model comes ONLY
%    from the post-encTime settling window (drive off, pure
%    self-excitation vs inhibition) -- during the drive-on phase
%    (duration <= EncTime), longer duration recruits MORE units, not
%    fewer, because nothing is yet suppressing the slowly-self-excited
%    tail units. So NewCommit(item 1) will likely still INCREASE from
%    Duration=0.1 to 0.4, then drop sharply at 0.8 -- a cliff, not the
%    smooth decrease you're after. I could not find a parameterization
%    within {inputDrive...asyW} that fixes this; it looks structural
%    to how EncTime interacts with Duration, not something commit/
%    inhib/asyW alone can smooth out. Try shrinking EncTime (e.g. to
%    ~0.15) only in combination with retesting commit/inhib -- in my
%    testing this tends to trade the problem for an even worse one
%    (item 1 can permanently block items 2-6 at every duration, not
%    just the longest one), so I'm not recommending it outright.
% ===================================================================

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
kappaB = 15; 
kappaCat = 10;
strengthSD = 0.1;
dnoise = 0.01;
nCat = 8;
nb = 90; 
encRate = 10;
primacy = 0.01;  % rate of exponential decline of encoding strength (smaller -> flatter)
commit = 0.70;   % <-- retuned from 0.85; see note above. Try 0.6-0.75 to explore the trade-off.
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
    inhib = 0.35;  
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
NewCommit = NaN(length(Duration), setsize);   % <-- NEW: mean newly-committed units per item

for dIdx = 1:length(Duration)
    duration = Duration(dIdx);
    encTime = min(duration, EncTime);
    nSteps = round(duration./tstep);

    Error = zeros(nTrials, setsize);
    maxW = zeros(nTrials, nSteps);
    NewCommitTrial = zeros(nTrials, setsize);   % <-- NEW

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

            prevCommittedMask = Committed > 0;               % <-- NEW
            Committed(W > cThreshold) = commit;
            newMask = (W > cThreshold) & ~prevCommittedMask;  % <-- NEW
            NewCommitTrial(trial, item) = sum(newMask(:));    % <-- NEW

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
    NewCommit(dIdx, 1:setsize) = mean(NewCommitTrial);   % <-- NEW

end

PreFigure([], [], 2);
subplot(1,2,1);
plot(1:setsize, SPC');
PostFigure([0.5, setsize+0.5, 0, 90], 'Input Position', 'Error', [], vec2legend(Duration));
subplot(1,2,2);
plot(Duration, meanError);
PostFigure([0, max(Duration)+0.1, 0, 90], 'Duration', 'Error');

% --- NEW diagnostic figure: newly-committed W-units by serial position,
% one line per duration. Plain MATLAB calls used here (not
% PreFigure/PostFigure) since I don't have those function definitions
% to match their exact call signature -- feel free to restyle to match
% your other figures.
figure;
plot(1:setsize, NewCommit', '-o');
xlabel('Input Position');
ylabel('Newly committed W-units');
legend(vec2legend(Duration));
title('Newly committed units per item, by duration');
xlim([0.5, setsize+0.5]);

% MaxW = mean(maxW);
% Time = tstep*(1:nSteps);
% subplot(2,3,setsize);
% plot(Time, MaxW);
% PostFigure([0, duration, 0, max(0.1, max(MaxW))], 'Time (s)', 'max(W)', num2str(setsize) );
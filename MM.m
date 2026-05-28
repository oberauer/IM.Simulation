function output = MM(modelParms, Data, mode)
% mixture model to be plugged into FitGMM wrapper (by HY, modified by KO)

maxSZ = size(Data.D,2);

sd = 2*pi*modelParms(1)/360;  % standard deviation in the von Mises distribution, translated from degrees to radians
kappa = 1/(sd^2);             %kappa of von Mises distribution, in radians
beta = modelParms(2); % probability of random guessing
if length(modelParms) > 2, gamma = modelParms(3); else, gamma = 0; end  %probability of transposition
if length(modelParms) > 3, delta = modelParms(4); else, delta = 0; end  %probability of color-wheel attraction effect

positionInd = repmat(1:maxSZ, size(Data.D,1), 1);
nonExistingItems = positionInd > repmat(Data.setsize, 1, maxSZ);  %index into items that don't exist for a given setsize
Data.D(nonExistingItems) = -1;

cuedItem = Data.D == 0;  %cuedItem = logical index into column of colors to identify the cued (target) color
transItem = Data.D >= .5;  %logical indices into non-target colors

transItem = transItem ./ repmat(sum(transItem, 2), 1, maxSZ);  %divide by number of non-targets (strictly, speaking, not necessary, because TP is normalized below anyway)
transItem(isnan(transItem)) = 0;   %nans can occur from division by zero at setsize 1

Pchoose = VonMises(Data.Dcang, 0, kappa); %for each trial and each item j in that trial, probability of choosing color i, given item j is retrieved

% Probability distribution centered on cued item
CStrength = Pchoose .* repmat(cuedItem, [1 1 360]);
CP = sum(CStrength, 2) ./ repmat(sum(sum(CStrength, 2), 3), [1, 1, 360]);  %normalize
CP = squeeze(CP);

% Probability distribution(s) centered on non-target items
TStrength = Pchoose .* repmat(transItem, [1 1 360]);
TP = sum(TStrength, 2) ./ repmat(sum(sum(TStrength, 2), 3), [1, 1, 360]);  %normalize
TP(isnan(TP)) = 0;
TP = squeeze(TP);

% Probability distribution centered on color-wheel color closest to target
CWStrength = VonMises(Data.Dwang, 0, kappa);
CWP = CWStrength./repmat(sum(CWStrength,2), 1, 360); 

guessing = ones(size(CP)) * 1 / 360;

pred = CP .* (1-beta-gamma-delta) + gamma .* TP + guessing .* beta + CWP .* delta;
if mode == 1, output = pred; end
if mode == 2
    I = repmat(1:360, size(Data.D,1), 1);
    D = repmat(Data.response, 1, 360);
    Index = I == D;                %for each trial (rows), Index is 1 in the xth column, where x is the actually selected color in that trial, and zero otherwise
    likelihood = pred(Index);      %Data codes the responses (i.e., the selected color angle), serving as index into the prob-distribution
    if min(likelihood) < 0, logL = 7777777; else logL = -sum(2*log(likelihood)); end %safeguard because in MM2, likelihoods can go negative if Ci+Ai > 1
    if isnan(logL), logL = 888888888; end
    if beta+gamma > 1, logL = 999999999; end
    output = logL;
end

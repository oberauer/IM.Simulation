function meanDiff = ResponseDistrib(Array, Target, Response)
% Plots response distributions around correct color and color of
% non-targets

nBin = 45;
szBin = 360 / nBin;

ub = (-180+szBin):szBin:180;
lb = -180:szBin:(180-szBin);
binLabel = (-180+szBin/2):szBin:(180-szBin/2);

Setsize = sum(Array > 0, 2);
setsizes = sort(unique(Setsize));
if size(setsizes,1) > size(setsizes,2), setsizes = setsizes'; end  % make sure that resulting vector is horizontal

Diff = Response - Target;
Diff(Diff>180) = Diff(Diff>180) - 360;
Diff(Diff<=-180) = Diff(Diff<=-180) + 360;
meanDiff = zeros(1, max(unique(Setsize)));

figure
fIdx = 1;
for setsize = setsizes
    distribution(setsize).respDist = zeros(1,nBin);
    diff = Diff(Setsize==setsize);
    meanDiff(setsize) = mean(abs(diff));
    for tIndex = 1:length(diff)
        distribution(setsize).respDist((diff(tIndex) >= lb) & (diff(tIndex) < ub)) = distribution(setsize).respDist((diff(tIndex) >= lb) & (diff(tIndex) < ub)) + 1;
    end

    subplot(2, ceil(length(setsizes)/2), fIdx);
    plot(binLabel, distribution(setsize).respDist ./ sum(distribution(setsize).respDist), '-k')
    ylim([0 .4]);
    title(['Setsize = ', mat2str(setsize), '; Errors around Target']);
    fIdx = fIdx + 1;
end
set (gcf, 'Color','w');
set (gca, 'box', 'on');


figure
fIdx = 1;
for setsize = setdiff(setsizes, 1)
    distribution(setsize).respDist = zeros(1,nBin);
    response = Response(Setsize==setsize);
    target = Target(Setsize==setsize);
    array = Array(Setsize==setsize, :);

    for tIndex = 1:length(response)
        tmpDistribution = zeros(1,nBin);
        ntargets = setdiff(array(tIndex,1:setsize), target(tIndex));
        for i = 2:setsize
            diff = ntargets(i-1) - response(tIndex);
            if diff > 180, diff = diff - 360; end
            if diff <= -180, diff = diff + 360; end
            tmpDistribution((diff >= lb) & (diff < ub)) = tmpDistribution((diff >= lb) & (diff < ub)) + 1;
        end
        tmpDistribution = tmpDistribution ./ (setsize-1); %divide by number of non-target memory items that entered into successive diff calculcations and counts in distribution
        distribution(setsize).respDist = distribution(setsize).respDist + tmpDistribution;
    end
    subplot(2, ceil(length(setsizes)/2), fIdx);
    plot(binLabel, distribution(setsize).respDist ./ sum(distribution(setsize).respDist), '-k');
    ylim([0 .1]);
    title(['Setsize = ', mat2str(setsize), '; Deviations from Nontargets']);
    fIdx = fIdx + 1;
end
set (gcf, 'Color','w');
set (gca, 'box', 'on');

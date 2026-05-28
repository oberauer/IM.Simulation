function ResponseDistribC(Array, Target, Response, Conditions, CondLabels)
% Plots response distributions around correct color, color of non-targets,
% and color in the color wheel closest to the target, for fixed set size
% and several experimental conditions

if (nargin < 5), CondLabels = mat2str(unique(Conditions)); end

nBin = 45;
szBin = 360 / nBin;
ub = (-180+szBin):szBin:180;
lb = -180:szBin:(180-szBin);
binLabel = (-180+szBin/2):szBin:(180-szBin/2);

setsize = size(Array, 2) - 1; 

Diff = wrap(Response - Target, 180);

% response distributions centered on the target
figure;
for cond = 1:max(Conditions)
    distribution(cond).respDist = zeros(1,nBin);
    diff = Diff(Conditions==cond);
    for tIndex = 1:length(diff)
        distribution(cond).respDist((diff(tIndex) >= lb) & (diff(tIndex) < ub)) = distribution(cond).respDist((diff(tIndex) >= lb) & (diff(tIndex) < ub)) + 1;
    end
    
    subplot(2, ceil(max(Conditions)/2), cond);
    plot(binLabel, distribution(cond).respDist ./ sum(distribution(cond).respDist), '-k')
    ylim([0 .4]);
    title('Errors; ', [CondLabels{cond}]);
end
set (gcf, 'Color','w');
set (gca, 'box', 'on');

% response distributions centered on the non-targetrs
figure;
for cond = 1:max(Conditions)
    distribution(cond).respDist = zeros(1,nBin);
    response = Response(Conditions==cond);
    target = Target(Conditions==cond);
    array = Array(Conditions==cond, :);
    
    for tIndex = 1:length(response)
        tmpDistribution = zeros(1,nBin);
        ntargets = setdiff(array(tIndex,1:setsize), target(tIndex));
        for i = 2:setsize
            diff = wrap(ntargets(i-1) - response(tIndex), 180);
            tmpDistribution((diff >= lb) & (diff < ub)) = tmpDistribution((diff >= lb) & (diff < ub)) + 1;
        end
        tmpDistribution = tmpDistribution ./ (setsize-1); %divide by number of non-target memory items that entered into successive diff calculcations and counts in distribution
        distribution(cond).respDist = distribution(cond).respDist + tmpDistribution;
    end
    subplot(2, ceil(max(Conditions)/2), cond);
    plot(binLabel, distribution(cond).respDist ./ sum(distribution(cond).respDist), '-k');
    ylim([0 .1]);
    title(['Dev. from NT; ', CondLabels{cond}]);
end
set (gcf, 'Color','w');
set (gca, 'box', 'on');


CWdiff = wrap(Response' - Array(:, setsize+1), 180); %Color in color wheel closest to the target location is in Array, column setsize+1

figure
for cond = 1:max(Conditions)
    distribution(cond).respDist = zeros(1,nBin);
    cwdiff = CWdiff(Conditions==cond);
    for tIndex = 1:length(cwdiff)
        distribution(cond).respDist((cwdiff(tIndex) >= lb) & (cwdiff(tIndex) < ub)) = distribution(cond).respDist((cwdiff(tIndex) >= lb) & (cwdiff(tIndex) < ub)) + 1;
    end
    
    subplot(2, ceil(max(Conditions)/2), cond);
    plot(binLabel, distribution(cond).respDist ./ sum(distribution(cond).respDist), '-k')
    ylim([0 .1]);
    title(['Dev. from CW; ', CondLabels{cond}]);
end
set (gcf, 'Color','w');
set (gca, 'box', 'on');

function [meanCTF] = ApplyIEM(W, EEG, Angle, channelCenters, itemIdx)
% trains and applies the Inverted Encoding Model, computes the channel
% tuning function from the EEG signal
% Angle must be given in degrees!

global E

% if no index matrix is given, use all stimuli
if nargin < 5
    itemIdx = repmat(1:size(Angle,2), E.ntrials, 1);
end

nChannels = size(W, 1);
chDeg = round(channelCenters);  

% extend scale to the left and the right to wrap around for interpolation
extendLeft = (max(chDeg):360) - 360;  % degrees to be added on the left of 1 degree: space from largest channel center to 360
extendRight = 360 + (1:min(chDeg));   % degrees to be added on the right of 360 degree: space from 1 to smallest channel center
xDeg = [extendLeft, 1:360, extendRight]; % x values for interpolation in degres
chDegExt = [extendLeft(1), chDeg, extendRight(end)];  % channel centers in degrees, extended by the right-most on the left, and the left-most on the right, for wrap-around

C0 = 180; % central channel, on which the CTFs will be centered

% Use inverted W to reconstruct channel responses for the second half of trials
ChannelResponse = (inv(W*W') * W * EEG').';

% now center each channelResponse on each of the array stimuli considered in each trial
nItems = size(itemIdx,2); % number of array stimuli considered for each trial
ChannelResponseC = zeros(size(EEG,1), nItems, nChannels);

for trial = 1:size(EEG,1)
    for idx = 1:nItems
        chResponse = ChannelResponse(trial,:); 
        channelResponse = [chResponse(end), chResponse, chResponse(1)]; % wrap-around for interpolation: extend to the left with the rightmost value, and to the right with the left-most value
        chResponse360 = interp1(chDegExt, channelResponse, xDeg); % interpolate linearly the 360 values on the Xq scale
        channelResponse360 = chResponse360(xDeg>0 & xDeg<361); % cut out the part that is not extended over the edges
        centeredResponse360 = circshift(channelResponse360, C0-round(Angle(trial, itemIdx(trial, idx))), 2); 
        ChannelResponseC(trial,idx,:) = centeredResponse360(chDeg);  % take the values of the centered response at the channel centers
    end
end
if (nItems == 1)
    meanCTF = squeeze(mean(ChannelResponseC, 1))';
else
    meanCTF = squeeze(mean(ChannelResponseC, 1));
end
end


function [meanCTF] = IEM(StimMask, basisSet, EEG, Pangle, itemIdx)
% trains and applies the Inverted Encoding Model, computes the channel
% tuning function from the EEG signal

global E

% if no index matrix is given, use all stimuli
if nargin < 5
    itemIdx = repmat(1:size(Pangle,2), E.ntrials, 1);
end

nChannels = size(basisSet, 2);
C0 = ceil(nChannels/2);  % the number of the central channel, on which the CTFs will be centered

% Train the IEM on the first half of trials
trainIdx = 1:(E.ntrials/2);
testIdx = (E.ntrials/2+1):E.ntrials;

% Design matrix: For each trial, create the predicted response of
% the channels, given the stimulus (stimuli) actually presented
PredChannelResponse = StimMask(trainIdx,:) * basisSet;

% Find the weight matrix projecting hypothetical channel responses to electrodes
W = PredChannelResponse\EEG(trainIdx, :);  % EEG = PredChannelResponse * W, solve for W

% Use inverted W to reconstruct channel responses for the second half of trials
ChannelResponse = (inv(W*W') * W * EEG(testIdx,:)').';

% now center each channelResponse on each of the array stimuli considered in each trial
nItems = size(itemIdx,2); % number of array stimuli considered for each trial
ChannelResponseC = zeros(E.ntrials/2, nItems, nChannels);

for trial = 1:length(testIdx)
    for idx = 1:nItems
        ChannelResponseC(trial,idx,:) =  circshift(ChannelResponse(trial,:), C0-Pangle(testIdx(trial), itemIdx(testIdx(trial),idx)), 2);
    end
end
if (nItems == 1)
    meanCTF = squeeze(mean(ChannelResponseC, 1))';
else
    meanCTF = squeeze(mean(ChannelResponseC, 1));
end
end


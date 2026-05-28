function W = TrainIEM(StimMask, basisSet, EEG)
% trains the Inverted Encoding Model on a training set of EEG patterns

% Design matrix: For each trial, create the predicted response of
% the channels, given the stimulus (stimuli) actually presented
PredChannelResponse = StimMask * basisSet;

% Find the weight matrix projecting hypothetical channel responses to electrodes
W = PredChannelResponse\EEG;  % EEG = PredChannelResponse * W, solve for W
% 
% EEGpred = PredChannelResponse * W;
% disp(EEGpred);

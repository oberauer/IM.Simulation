function retrieved = ccorr(cue, memory, n)
% ccorr(retrieval cue, circular-convolution memory, number of elements)
% circular correlation - the inverse of circular convolution
invcue = zeros(size(cue));
invcue(1) = cue(1);
invcue(2:end) = fliplr(cue(2:end));  % involution
retrieved = cconv(memory, invcue, n);    % retrieve 



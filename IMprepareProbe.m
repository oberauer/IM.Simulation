function [probestim, probeIdx] = IMprepareProbe(F, setsize)
% prepares probes for recognition test

global E
global C

if E.test == 2   % old-new recognition
    if E.ptype == 1, probeIdx = F(1,1); end   % select the probe color for positive probes (equal to the first array item, which is always the target)
    if E.ptype == 2, probeIdx = F(1,setsize+1); end  % select the color for new probes (randomly from all possible features in the candidate set)
    if E.ptype == 3, probeIdx = F(1,randperm(setsize-1, 1)+1); end  % select the color of intrusion probes (randomly from any non-target, that is, from array items 2 to N
    probestim = C.stim(probeIdx,:);  % define the probe stimulus
end
if E.test == 3  % n-AFC
    probeIdx = zeros(1, length(E.respAlt));
    probestim = zeros(length(E.respAlt), C.nc); 
    probeIdx(1) = F(1, E.respAlt(1)); % usually this is the positive (correct) alternative
    probestim(1,:) = C.stim(probeIdx(1),:);
    for j = 1:sum(E.respAlt==2)
        probeIdx(j+1) = F(1, setsize+j); % extra-set probes
        probestim(j+1,:) = C.stim(probeIdx(j+1),:);
    end
    for k = 1:sum(E.respAlt==3)   
        probeIdx(j+k+1) = F(1, k+1);  % intrusion probes
        probestim(j+k+1,:) = C.stim(probeIdx(j+k+1),:);
    end
end
if E.test == 4   % change localization
    probestim = zeros(E.rss, C.nc);
    probeIdx = zeros(1, E.rss);  % here, probeIdx is not the feature value (F) but the index INTO F and L, because it needs to also code the location of the probes!
    lures = 1 + randperm(setsize-1); % shuffle the non-targets
    if E.ptype == 2, probeIdx(1) = setsize+1; probestim(1,:) = C.stim(F(1, probeIdx(1)), :); end % replace the target with the first new item in the F list
    if E.ptype == 3, probeIdx(1) = lures(end); probestim(1,:) = C.stim(F(1, probeIdx(1)), :); end % replace the target by the last in the shuffled item list
    for j = 2:E.rss
        probeIdx(j) = lures(j);
        probestim(j,:) = C.stim(F(1, probeIdx(j)),:);  % the remainder of the response set remains unchanged
    end
end

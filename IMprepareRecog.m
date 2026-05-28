function IMprepareRecog(expectedChange)
% prepares decision criterion for recognition

if nargin < 1, expectedChange = 0; end

global E
global C
global P

% prepare recognition decision criterion
likSame = VonMisesN(C.x, pi, P.kappacrit*P.kappa_feat);  % estimated likelihood of "same" trials using the meta-cognitive estimate of feature precision

if expectedChange == 0  % all stimuli equally likely as change probes
    featurestep = floor(C.nc/C.nstim);
    Lures = deg2rad(setdiff(0:featurestep:359, 180)); 
    likChange = mean(VonMisesN(C.x, Lures, P.kappacrit*P.kappa_feat));
    E.SameRange = likSame > likChange; 
else
    likChange1 = VonMisesN(C.x, pi+deg2rad(expectedChange), P.kappacrit*P.kappa_feat); % expect large changes: mean 90 degrees
    likChange2 = VonMisesN(C.x, pi-deg2rad(expectedChange), P.kappacrit*P.kappa_feat); % expect large changes: mean 90 degrees
    E.SameRange = (likSame > likChange1) .* (likSame > likChange2);  % range of retrieved feature values that are similar enough to the probe to say "same"
end



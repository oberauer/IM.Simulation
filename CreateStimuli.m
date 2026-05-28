function [maskFeatures] = CreateStimuli(nMaskStim)
% Function for creating stimuli

global C
global P
global E

if (nargin < 1), nMaskStim = min(C.nstim, 100); end  % default: mask with a large number of stimuli

% meanCtime = -(log(1-P.cTau)./P.cRate);
% C.meanNsteps = meanCtime/C.tstep;

C.stim = zeros(C.nstim, C.nc);
C.feature = zeros(C.nstim, 1);
featurestep = floor(C.nc/C.nstim);
if E.material == 3, featurestep = 1; end
for stimIdx = 1:C.nstim
    C.feature(stimIdx,:) = stimIdx*featurestep;
    C.stim(stimIdx,:) = C.ContentFun(C.x, deg2rad(stimIdx*featurestep), P.kappaf_feat);
end
if E.material == 3  % oriented bars: add a 2nd peak 180 degrees of the original peak
    for stimIdx = (C.nstim+1):(2*C.nstim)
        C.feature(stimIdx,:) = stimIdx*featurestep;
        C.stim(stimIdx-C.nstim,:) = C.stim(stimIdx-C.nstim,:) + C.ContentFun(C.x, deg2rad(stimIdx*featurestep), P.kappaf_feat);
    end
end
C.stimnoise = zeros(1, C.nCat);

maskFeatures = randperm(C.nstim, nMaskStim);
%C.maskStim = mean(C.stim(randperm(C.nstim, nMaskStim), :));  % take a random subset of n stimuli (e.g., 8 in Ricker & Hardmann, 2017), and average them (~addition with divisive normalization for limited attention)
C.maskStim = sum(C.stim(maskFeatures, :));  % take a random subset of n stimuli (e.g., 8 in Ricker & Hardmann, 2017), and average them (~addition with divisive normalization for limited attention)



% Creation of the location vectors (population codes )
C.Location = linspace(0, C.nc, C.nloc+1);  % nloc equally spaced locations (location nloc+1 is nc, so virtually identical to 1 on the circle, and is ignored)
C.Location = C.Location(1:C.nloc) + (C.Location(2)-C.Location(1))/2;
C.location = zeros(C.nloc, C.nc);
for loc = 1:C.nloc
    C.location(loc, :) = C.ContextFun(C.x, deg2rad(C.Location(loc)), P.kappaf_ctx);
end
C.locationnoise = P.a * ones(1, C.nLocCat);  % context background noise

% Compute default CW matrix (to be rotated for each trial)
C.euclidDist = sqrt((P.rad1 - cos(C.x)).^2 + sin(C.x).^2);% euclidean distance of each color-wheel color C.x from location 0
targetLocAngle = 0;
CWcolor = 0;
CW = zeros(C.nc);
CWweight = C.ContextFun(C.euclidDist, 0, P.kappaf_ctx); % weight of encoding of each color-wheel color as a function of its distance from location-angle 0 degrees
for loc = 1:C.nc
    deltaLoc = wrap(loc - targetLocAngle, 180);  % deviation of loc from target location
    colorAtLoc = CWcolor + deltaLoc; % as we shift the location away from the target location, the color-wheel color closest to it shifts in the same direction
    colorAtLoc(colorAtLoc < 1) = colorAtLoc + 360;
    colorAtLoc(colorAtLoc > 360) = colorAtLoc - 360;
    CW = CW + circshift(CWweight, loc, 2)' * C.ContentFun(C.x, deg2rad(colorAtLoc), P.kappaf_feat);
    % the color in this loc is encoded by associating it to CWweight, which is the distribution of spatial attention to all locations in the color wheel, given that attention is focused on the target location
    %CW = CW + CWweight(loc) * (C.ContextFun(C.x, deg2rad(loc), P.kappaf_ctx)' * C.ContentFun(C.x, deg2rad(colorAtLoc), P.kappaf_feat));
end
C.CW = CW;  % strength of color wheel is reduced by lateral inhibition between the colors in close proximity

% distance matrix in FX
C.distFX = zeros(360);
for loc = 1:360
    C.distFX(loc,:) = abs(wrap(((1:360) - loc), 180)); 
end


end



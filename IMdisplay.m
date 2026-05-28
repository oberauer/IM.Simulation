function [Map, CWcolor, eraseFX] = IMdisplay(Map, L, probestim, probeIdx)
% adds the rotated color wheel, or the probe stimulus, to the display

global C
global E
global P

if nargin < 4, probeIdx = 0; end

Testdisplay = zeros(C.nc);
if E.test == 1   % continuous reproduction
    CWcolor = randperm(C.nstim,1);  % color of the color wheel closest to the target location at test
    if E.wheel == 1
        [Map, eraseFX] = UpdateFX(Map, 1);  % reset of target feature map because a new stimulus belonging to that feature map (ie, the color wheel) is attended
        targetLocAngle = find(C.location(L(1),:) == max(C.location(L(1),:)));
        Testdisplay = circshift(C.CW, [targetLocAngle, CWcolor]);
    else
        eraseFX = 0;
    end
end 
if E.test == 2  % change detection (yes/no)
    [Map, eraseFX] = UpdateFX(Map, 1);  % reset of target feature map because a new stimulus belonging to that feature map (probe) is attended
    Testdisplay = Testdisplay + C.location(L(1),:)' * probestim; % encode the probe at the target location
    CWcolor = 0; 
end
if E.test == 3  % change localization
    [Map, eraseFX] = UpdateFX(Map, 1);  % reset of target feature map because a new stimulus belonging to that feature map (probe) is attended
    Testdisplay = Testdisplay + C.location(L(1),:)' * probestim(1,:);  % target is L(1)
    for probe = 2:E.rss
        Testdisplay = Testdisplay + C.location(L(probeIdx(probe)),:)' * probestim(probe,:);  
    end
    CWcolor = 0; 
end

Map(1).FX = Map(1).FX + eraseFX * P.filter(E.test) * Testdisplay; %add test display with color wheel or probe to the target-feature map


end
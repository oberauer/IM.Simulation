function [Data] = Dataprocessing(Probedpos, pangle, cangle, target, response, setsize, Colors)
% Pre-processing of data from VisWM BindingMatrix simulations for fitting
% with GMM models - adapted from data pre-processing in the FitGMM code

sz = size(cangle,2)-1; 
Data.setsize = setsize; 
Data.response = response;
for t = 1:size(Probedpos,1)
    Data.probed(t) = find(Probedpos(t) == pangle(t,:), 1);
end
Data.Tang = target;
Data.Cang = cangle(:, 1:sz); 

probed = repmat(Probedpos, 1, size(pangle,2));
D1 = abs(pangle-probed);         %distance of probed location to all 6 locations
D2 = abs((pangle-13)-probed);    %distance of probed location to all 6 locations, setting 13 to 0 etc, to capture distance across the 13->1 boundary
D3 = abs(pangle - (probed-13));  %ditto for distance between large Probed and smaller pangle)
Data.D = min(min(D1, D2), D3);   %keep smallest distance
Data.D = 2*pi*Data.D./13;        %re-normalize to radians

Cangle = repmat(cangle(:,1:sz), [1, 1, 360]);
Dc1 = abs(Cangle-Colors);
Dc2 = abs((Cangle-360)-Colors);
Dc3 = abs(Cangle-(Colors-360));
Dc = min(min(Dc1, Dc2), Dc3);   %distances of each item in each trial from the 360 colors in color wheel: matrix of N x setsize x 360
Data.Dcang = 2*pi*Dc/360;       %convert into scale from 0 to 2pi (for 1:360), so that distances range from 0 to pi.

Tangle = repmat(target, [1, 360]);
Dt1 = abs(Tangle - squeeze(Colors(:,1,:)));
Dt2 = abs((Tangle-360) - squeeze(Colors(:,1,:)));
Dt3 = abs(Tangle - (squeeze(Colors(:,1,:))-360));
Dt = min(min(Dt1, Dt2), Dt3);
Data.Dtang = 2*pi*Dt/360;       %distance of target item in each trial from the 360 colors in color wheel (transformed to radians)

CWangle = repmat(cangle(:,sz+1), [1, 360]);   %color-wheel color closest to target is coded in cangle column "setsize+1"
Dw1 = abs(CWangle - squeeze(Colors(:,1,:)));
Dw2 = abs((CWangle-360) - squeeze(Colors(:,1,:)));
Dw3 = abs(CWangle - (squeeze(Colors(:,1,:))-360));
Dw = min(min(Dw1, Dw2), Dw3);
Data.Dwang = 2*pi*Dw/360;       %distance of color-wheel color closest to target in each trial from the 360 colors in color wheel (transformed to radians)


function [Map, W, GateClosed, gWeight, Focus, Afocus, content, context, SpatAttn] = IMdistraction(Map, W, GateClosed, gWeight, L, F, nDistr, tDistr)
% encodes a memory set, simultaneously or sequentially

global P
global E
global C

SpatAttn = zeros(C.nc, 1);
content = zeros(E.nfeat, C.nCat);
Afocus = zeros(E.nfeat, C.nc);
for distr = 1:nDistr
    Focus = distr;
    Map = UpdateFX(Map);  % complete reset of feature map because a new stimulus is attended
    for ff = 1:E.nfeat
        Map(ff).FX = Map(ff).FX + C.location(L(Focus),:)' * (C.stim(F(ff,Focus),:));
    end
    % focus attention to the relevant stimulus location in the feature map(s) FX
    AfocusLoc = C.location(L(Focus),:);     % update location attended to in the feature maps
    for ff = 1:E.nfeat, Afocus(ff,:) = AfocusLoc./sum(AfocusLoc) * Map(ff).FX; end % use location as (spatial) attentional filter to pull out the target feature from its feature map
    if E.context == 1, context = C.locationnoise + AfocusLoc * C.MappingC; end  %
    if E.context == 2, context = C.stimnoise + (AfocusLoc./sum(AfocusLoc) * Map(2).FX) * C.Mapping; end
    for ff = 1:E.nfeat, content(ff,:) = C.stimnoise + Afocus(ff,:) * C.Mapping; end
    Map(1).FX = Map(1).FX + C.tstep * (-Map(1).FX + Map(1).FX .* repmat(SpatAttn, 1, C.nc)); % spatial attention modulates feature maps
    [W, GateClosed, gWeight, ~, ~] = IMencodeStim(W, context, content, GateClosed, gWeight, P.cRate, tDistr); % consolidate distractor in WM-binding module
end





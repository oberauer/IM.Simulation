function [Map, eraseFX] = UpdateFX(Map, dimension)

global P
global C

if nargin < 2, dimension = 1:C.nfeatures; end
% by default, update all dimensions in the map

eraseFX = rand > P.eraseFX; 
for fidx = intersect(dimension, (1:C.nfeatures))
    Map(fidx).FX = (1-eraseFX)*Map(fidx).FX;  % if test display is attended, it erases FX (and encodes new stimulus)  
end

end
function [basisSet, eW, eW2, nElectrodes, channelCenters] = CreateIEM(nElectrodes, nChannels, eRegular, eGrad,  plotting)
%Set up an Inverted Encoding Model

global C
global P

if nargin < 5, plotting = 0; end

% create a square-head EEG electrode distribution

eegx = round(sqrt(nElectrodes));  % number of x coordinates
eegy = ceil(nElectrodes/eegx);    % number of y coordinates
nElectrodes = eegy*eegx;  % ... then update to the number of electrodes in a rectangular matrix
[eX, eY] = meshgrid(1:eegx, 1:eegy);  % actual coordinates of all electrodes
eX = eX - mean(mean(eX));  % center on mean, because x-y coordinates of locations are also given relative to screen center as [0,0]
eY = eY - mean(mean(eY)); 

eW = 0.3 + randn(C.nc, nElectrodes); % projection matrix from 360 locations on a virtual circle to the electrodes: initialize some randomness; add 0.3 for global shift -> CDA (negative = up!)
eW2 = 0.3 + randn(P.nb, nElectrodes); % mapping from binding layer to electrodes 
xLoc = cosd(1:360)*(0.3*eegx);  % radius set to 0.3 of the extension of the electrode matrix
yLoc = cosd(1:360)*(0.3*eegy);
for loc = 1:C.nc
    for elec = 1:nElectrodes
        D = sqrt((xLoc(loc)-eX(elec))^2 + (yLoc(loc)-eY(elec))^2);  % distance of electrode from the location
        eW(loc, elec) = eW(loc, elec) + eRegular * exp(-eGrad*D);   % add some regular spatial mapping from screen location to electrode location
    end
end

if plotting
    % plot the EEG profile (over the 50 electrodes) for each of the locations
    PreFigure;
    for L = 1:length(C.Location)
        subplot(4,4,L);
        loc = round(C.Location(L));
        scalpProfile = reshape(eW(loc,:), eegy, eegx);
        imagesc(scalpProfile);
    end
end

channelCenters = 360*(0.5:(nChannels+0.49))./nChannels; % channel centers on the degree scale
make_basis_function = @(xx,mu) (cosd((xx-mu)/2)).^(nChannels-mod(nChannels,2));  % function for the shape of the channel tuning functions
basisSet = zeros(C.nc, nChannels);   % projection of locations onto the hypothetical channels
for cc = 1:nChannels
    basisSet(:,cc) = make_basis_function(1:C.nc,channelCenters(cc));
end

end


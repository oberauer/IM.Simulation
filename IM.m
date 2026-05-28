% Compute predictions for IM for Color Wheel experiment 1:
% This is the closed-form version of the IM published in Oberauer & Lin
% (2017)

function output = IM(modelParms, Data, mode)
C = 1;
B = modelParms(1);
A = modelParms(2);
s = modelParms(3);
kappa = modelParms(4);
kappafocus = modelParms(5);
Creduction = modelParms(6);

cueing = exp(-s*Data.D);      %strength of cue-based retrieval for each memory object, as a function of spatial distanc
positionInd = repmat(1:max(Data.setsize), length(Data.setsize), 1);
nonExistingItems = positionInd > repmat(Data.setsize, 1, max(Data.setsize));  %index into items that don't exist for a given setsize
cueing(nonExistingItems) = 0; %zero strength for objects that don't exist

Pchoose = VonMisesN(Data.Dcang, 0, kappa, 3); %for each trial and each item j in that trial, probability of choosing color i, given item j is retrieved (and is not focused)
Pchoose(repmat(nonExistingItems, [1 1 360])) = 0; % this is necessery for A component
Pchoosefocus = squeeze(VonMisesN(Data.Dcang(:,1,:), 0, kappafocus, 3)); %for each trial, probability of choosing color i, given the target is focused (high-precision representation exists only for one item, does not spill over to neighbors!)

N = size(Data.setsize,1);
Evidence = zeros(N, 360);
for trial = 1:N
    Eb = (1/360)*Data.setsize(trial);
    Ea = squeeze(sum(Pchoose(trial,:,:),2))';
    Ec = squeeze(sum(bsxfun(@times, cueing(trial,:), Pchoose(trial,:,:))))';
    if Data.preretro==1  % pre-cue, or no cueing (in which case Data.cueing is always = 1)
        switch Data.cueing(trial)
            case 1, pFocus = 1/Data.setsize(trial);   %neutral
            case 2, pFocus = 1;                         %valid
            case 3, pFocus = 0;                         %invalid
        end
        if Data.setsize==1, pFocus = 1; end             %even in invalid-cue case
        Evidence(trial,:) = (1-pFocus(trial)).*(B*Eb + A*Ea + C*Ec) + pFocus(trial).*(Creduction*B*Eb + Creduction*A*Ea + C*(Ec + Pchoosefocus(trial,:)));
    end
    if Data.preretro==2  % post-cue
        switch Data.cueing(trial)
            case 1, pfenc = 1; pfnenc = 0;   %neutral: P(focus at retrieval | focus at encoding) = 1, P(focus at retrieval | not focus at enc) = 0
            case 2, pfenc = 1; pfnenc = 1;   %valid
            case 3, pfenc = 0; pfnenc = 0;   %invalid
        end
        Pfe = 1./Data.setsize(trial); %probability of focusing the target at encoding
        pfocused = Pfe.*pfenc + (1-Pfe).*pfnenc; %probability of focusing the target at retrieval
        Evidence(trial,:) = (1-pfocused).*(Eb + A*Ea + C*Ec) ...
            + Pfe*pfenc * (Creduction*Eb + Creduction*A*Ea + C*(Ec + Pchoosefocus(trial,:)) ) ...
            + (1-Pfe)*pfnenc * (Creduction*Eb + Creduction*A*Ea + C*Ec);
    end
end
pred = bsxfun(@rdivide, Evidence, sum(Evidence, 2));  %normalize by dividing by the sum across all 360 colors
if mode == 1, output = pred; end
if mode == 2
    I = repmat(1:360, size(Data.D,1), 1);
    D = repmat(Data.response, 1, 360);
    Index = I == D;                %for each trial (rows), Index is 1 in the xth column, where x is the actually selected color in that trial, and zero otherwise
    likelihood = pred(Index);      %Data codes the responses (i.e., the selected color angle), serving as index into the prob-distribution
    if min(likelihood) < 0, logL = 7777777; else logL = -sum(2*log(likelihood)); end %safeguard because in MM2, likelihoods can go negative if Ci+Ai > 1
    if isnan(logL), logL = 888888888; end
    output = logL;
end


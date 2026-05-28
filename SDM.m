% compute predictions for the full 4-parameter SDM

function output = SDM(modelParms, Data, mode)

%x = pi*(-179:180)./180;  % all possible features, in radians
x = pi*(1:360)./180;

C = modelParms(1);
kappa  = modelParms(2);
A = modelParms(3);
s = modelParms(4);

nTrials = length(Data.response);
Pred = zeros(nTrials, length(x));
%Pred2 = zeros(nTrials, length(x));
%L1 = zeros(nTrials, 1);
%L2 = zeros(nTrials, 1);

Signal = VonMises(Data.Dcang, 0, kappa); %for each trial and each item j in that trial, activation of color i through the cueing of item j
% Dcang codes, for each trial and each item in that trial, the distance of
% that item's feature from all 360 colors. Hence, Dcang = 0 at the item's
% feature value. The von

Cueing = C * exp(-s*Data.D);

for trial = 1:nTrials
    Act = (Cueing(trial,:) + A) * squeeze(Signal(trial,:,:));
    Pred(trial,:) = exp(Act)./sum(exp(Act));
    % if isnan(Pred(trial, Data.response(trial)))
    %     L1(trial) = eps;
    % else
    %     L1(trial) = Pred(trial, Data.response(trial));  % we need the response, not the error!
    % end
end

I = repmat(1:360, size(Data.D,1), 1);
D = repmat(Data.response, 1, 360);
Index = I == D;                %for each trial (rows), Index is 1 in the xth column, where x is the actually selected color in that trial, and zero otherwise
L = Pred(Index);      %Data codes the responses (i.e., the selected color angle), serving as index into the prob-distribution

% for trial = 1:nTrials
%     cueing = C * exp(-s*Data.D(trial,:));
%     % Data.Distance(trial,:) is the vector of spatial distances between
%     % each feature and the target (for the target, the distance is 0, of
%     % course)
%     Act = zeros(1, length(x));
%     for item = 1:Data.setsize(trial)
%         Act = Act + (cueing(item) + A) * VonMises(x, Data.feature(trial, item), kappa);
%         % Data.feature(trial, :) is the vector of features (in radians)
%         % that were presented in the current trial
%     end
%     Pred2(trial,:) = exp(Act)./sum(exp(Act));
%     if isnan(Pred2(trial, Data.response(trial)))
%         L2(trial) = eps;
%     else
%         L2(trial) = Pred2(trial, Data.response(trial));  % we need the response, not the error!
%     end
% end

%disp([sum(L), sum(L1), sum(L2)]);

if mode == 2
    Deviance = -2*sum(log(L));
    output = Deviance;
else
    output = Pred;
end





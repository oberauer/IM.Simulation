function ParX = CreateIndDiff
% Creates individual differences of parameter values

global P
global E
global C

if ~isempty(C.indVar)
    for ii = 1:length(C.indVar)
        eval(['meanIndVar(ii) = P.', C.indVar{ii}, ';']);
    end
end
if E.ngroups > 1
    for ii = 1:length(C.indVar)
        pGroup(:,ii) = exp( ((1:E.ngroups) - mean(1:E.ngroups)) .* C.SDfactor(ii) );  % proportion of mean parameter value to take as mean in each group
    end
else
    pGroup = ones(1,length(C.indVar));
end

ParX = zeros(E.nsubj*E.ngroups, length(C.indVar));

for group = 1:E.ngroups
    for ii = 1:length(C.indVar)
        if meanIndVar(ii) == C.nullVar(ii)
            ParX(((group-1)*E.nsubj+1):group*E.nsubj, ii) = meanIndVar(ii);  % don't vary across subjects if the mean parameter value is the null value of that parameter!
        else
            if C.logistVar(ii)==1
                maxVar = logit(C.maxIndVar(ii)); meanVar = logit(meanIndVar(ii));
                logitPar = min( maxVar, randn(1, E.nsubj)*C.SDfactor(ii)*meanVar + meanVar );
                ParX(((group-1)*E.nsubj+1):group*E.nsubj, ii) = logist(logitPar);
            else
                minIndVar = meanIndVar(ii) - 0.99*abs(meanIndVar(ii));
                ParX(((group-1)*E.nsubj+1):group*E.nsubj, ii) = min(C.maxIndVar(ii), max(minIndVar, meanIndVar(ii) * pGroup(group, ii) + randn(1, E.nsubj)*C.SDfactor(ii)*abs(meanIndVar(ii)) ));
            end
        end
    end
end
% make sure that nb is an integer
nbIdx = find(strcmp(C.indVar, 'nb'));
if ~isempty(nbIdx)
    ParX(:, nbIdx) = round(ParX(:, nbIdx));
end

E.nsubj = E.nsubj * E.ngroups;  % update number of subjects to actualy simulate

end


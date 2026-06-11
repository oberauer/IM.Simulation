%%% Toy model for investigating how CDA amplitude (from GateClosed)
%%% develops with set size

clear all
%close all

nb = 100;
delta = 0.8;
pMax = 1;
pBase = 0.3;
rRate = 4;

releaseTime = 0.5;  % that's about when the amplitude is usually measured

nRuns = 100;
CDA = zeros(nRuns, 8);
SPC = NaN(8);

for setsize = 1:8
    BStrength = zeros(nRuns, setsize);
    for run = 1:nRuns
        GateClosed = zeros(1, nb);
        BindingsMaintained = zeros(setsize, nb);
        for item = 1:setsize
            freeBP = rand(1, nb) < (1-delta);
            GateClosed(freeBP) = 0;
            BindingsMaintained(1:(item-1), freeBP) = 0; 
            initRecruited = rand(1, nb) < pMax;
            Binding = randn(1, nb)./sqrt(nb);   % normalize
            Binding = Binding .* (1-GateClosed) .* initRecruited;
            [Bsorted, sortIdx] = sort(abs(Binding), 'descend');
            ncommit = round( sum((1-GateClosed) .* initRecruited) * (pBase + (1-pBase) * exp(-rRate*releaseTime) ) );  % number of binding units that will remain committed -> gate will be closed
            committedNew = sortIdx(1:ncommit);
            GateClosed(committedNew) = 1;  % take the nbind free units with the highest absolute gating values and close the gate for them (= commit them)
            BindingsNow = zeros(1, nb);
            BindingsNow(committedNew) = Binding(committedNew);
            BindingsMaintained(item, :) = BindingsNow;
        end
        BStrength(run, :) = sum(abs(BindingsMaintained),2);
        CDA(run, setsize) = sum(GateClosed);

    end
    SPC(setsize, 1:setsize) = mean(BStrength);
end

PreFigure;
subplot(1,2,1);
plot(1:8, mean(CDA));
PostFigure([0.5, 8.5, 0, max(CDA(:))], 'Setsize', 'CDA');
subplot(1,2,2);
plot(1:8, SPC');
PostFigure([0.5, 8.5, 0, 8], 'Setsize', 'Binding Strength');

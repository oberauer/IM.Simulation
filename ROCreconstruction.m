%%% Test He, Kellen & Singmann's ROC reconstruction

m = 5;
R1 = zeros(1, m);  % R(5) computed by Equation 3
R2 = zeros(1, m);  % R(5) computed by Equation 2
PCorrect = [1.0000    0.6100    0.3750    0.3400    0.2350];
% proportion of correct choices for m = 1:5

% Equation 3
R1(1) = PCorrect(5);
R1(2) = 4*(PCorrect(4) - PCorrect(5));
R1(3) = 6*(PCorrect(3) - 2*PCorrect(4) + PCorrect(5));
R1(4) = 4*(PCorrect(2) - 3*PCorrect(3) + 3*PCorrect(4) - PCorrect(5));
R1(5) = PCorrect(1) - 4*PCorrect(2) + 6*PCorrect(3) - 4*PCorrect(4) + PCorrect(5);

% Equation 2
for i = 1:m
    Sum = 0;
    for j = (m-i+1):m
        Sum = Sum + nchoosek(i-1, j-(m-i+1)) * (-1).^(j-(m-i+1)) * PCorrect(j);
    end
    R2(i) = nchoosek(m-1, i-1) * Sum;
end

disp(R1);
disp(R2);

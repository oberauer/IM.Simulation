% Translation of R script BlockMarshakInequalities_fit_example.R by Claude

PCorrect = [0.610; 0.375; 0.340; 0.235]; % column vector
PCorrect = [0.6800    0.5050    0.4400    0.3300]; 

%% MM_5 matrix (14x6)
% R fills matrices column-by-column, same as MATLAB, so the structure is preserved.
% Columns: max | 2AFC | 3AFC | 4AFC | 5AFC | RHS
MM_5 = [
    2,  1,  0,  0,  0,  1/2;
    3,  0,  1,  0,  0,  1/3;
    4,  0,  0,  1,  0,  1/4;
    5,  0,  0,  0,  1,  1/5;
    2, -1,  0,  0,  0,   -1;
    3,  1, -1,  0,  0,    0;
    4,  0,  1, -1,  0,    0;
    5,  0,  0,  1, -1,    0;
    3, -2,  1,  0,  0,   -1;
    4,  1, -2,  1,  0,    0;
    5,  0,  1, -2,  1,    0;
    4, -3,  3, -1,  0,   -1;
    5,  1, -3,  3, -1,    0;
    5, -4,  6, -4,  1,   -1
];

%% Quadratic Programming
% R's QP.Solve(D, d, A, b) solves:  min  1/2 x'Dx - d'x
%                                   s.t. A'x >= b
%
% With A = -t(MM_5[,2:5]) and b = -MM_5[,6], the constraint becomes:
%   (-MM_5(:,2:5))' * x  >=  -MM_5(:,6)
%   =>  MM_5(:,2:5) * x  <=   MM_5(:,6)      (multiply both sides by -1)
%
% MATLAB's quadprog(H, f, A, b) solves:  min  1/2 x'Hx + f'x
%                                        s.t. A*x <= b
%
% Mapping:  H = eye(4),  f = -PCorrect,  A = MM_5(:,2:5),  b = MM_5(:,6)

% Claude's correction:
% Need: MM_5(:,2:5)*x >= MM_5(:,6)
% quadprog uses A*x <= b, so negate both sides:

H      = eye(4);
f      = -PCorrect;
A_ineq = -MM_5(:, 2:5);  % columns 2-5 (excludes "max" and "RHS")
b_ineq = -MM_5(:, 6);    % RHS column

x = quadprog(H, f, A_ineq, b_ineq);


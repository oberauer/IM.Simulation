function prob = VonMisesN(x, Mu, Kappa, xdim)
%VonMisesP.m (ko, modified from VonMises_quantiles)
% Compute (normalized) probabilities of histogram bins centered on values x for a vonMises distribution with mean mu and
% concentration Kappa (range of x: 0 to 2pi).  
% Arguments: x = a vector (or matrix) of equally spaced values between 0
% and 2*pi
% mu = mean or vector of means, Kappa = precision, 
% xdim = dimension of x along which the x-values are arranged (i.e.,
% dimension over which the sum for normalization runs)

if nargin < 4, xdim = 2; end

const = 2.*pi.*besseli(0,Kappa);
if length(Mu) > 1
    x = repmat(x, length(Mu), 1);
    mu = repmat(Mu', 1, length(x));  %in case mu is a vector, the output is a matrix with one row for each mu
else
    mu = Mu;
end
y = exp(Kappa.*cos(x-mu))./const;
prob = bsxfun(@rdivide, y, sum(y, xdim));  %normalize!

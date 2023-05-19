function W_mle = mle_fn(x,y,params)
% This function implements the MLE cost function for the 2D Gaussian.
%
% Inputs:
%   x: Nx2 matrix of x,y coordinates
%   y: Nx1 vector of counts
%   params: 5x1 vector of parameters
%
% Outputs:
%   W_mle: scalar value of the MLE cost function
%
%
% Created by Weihong Yeo, Northwestern University, 2022-12-29.
% Last modified by Weihong Yeo, Northwestern University, 2023-05-18.
%

a  = params(1);
b1 = params(2);
b2 = params(3);
c  = params(4);
d  = params(5);
mu_k = gaussian2dcirc(x(:,1),x(:,2),a,b1,b2,c,d+0.001);
n_k = y + 0.001;
W_mle = -sum(n_k .* log(mu_k) - mu_k - log(gamma(n_k)),'all');
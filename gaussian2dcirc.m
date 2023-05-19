function out = gaussian2dcirc(x,y,a,bx,by,c,d)
% This function returns a symmetric 2D Gaussian function with offset
% Implements the following function:
% 
% $$f(x,y) = a \exp\left[ -\frac{1}{2}\left( \left(\frac{x-b_x}{c_x} \right)^2 
% + \left(\frac{y-b_y}{c_y} \right)^2 \right)\right] + d$$
% 
% where $c_x = c_y = c$, such that the function is symmetric, and $x$ and $y$ 
% are assumed to be independent.
%
% Created by Weihong Yeo, Northwestern University, 2022-03-23.
% Last modified by Weihong Yeo, Northwestern University, 2022-03-23.
%

out = a.*exp(-0.5*((x-bx)/c).^2-0.5*((y-by)/c).^2)+d;

end
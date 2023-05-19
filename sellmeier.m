function eta = sellmeier(wl,glass)
% This function implements the Sellmeier equation for glass.
% calculates the refractive index of glass for given wavelengths.
% 
% Inputs:
%   wl: wavelength in units of μm
%   glass: glass type, either 'N-SF11' or 'N-LAF21', or provide the refractive
%          index in terms of an array:
%          glass = [B; C];
%          where B and C are 1x3 (row) vectors.
% Outputs:
%   eta: refractive index of glass
%
%
% Created by Weihong Yeo, Northwestern University, 2021-12-29.
% Last modified by Weihong Yeo, Northwestern University, 2023-05-18.
%

if strcmp(glass, 'N-SF11')
    B = [1.73759695, 0.313747346, 1.898781010];
    C = [0.01318870700, 0.0623068142, 155.2362900];
elseif strcmp(glass, 'N-LAF21')
    B = [1.87134529, 0.250783010, 1.220486390];
    C = [0.00933322280, 0.0345637762, 83.2404866];
elseif strcmp(glass, 'N-BK7')
    B = [1.03961212, 0.231792344, 1.01046945];
    C = [0.00600069867, 0.0200179144, 103.560653];
else
    B = glass(1,:);
    C = glass(2,:);
end
wl2 = wl .* wl;
eta = sqrt(1 + sum(B'*wl2 ./ (wl2 - C')));
end
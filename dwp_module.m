function pos = dwp_module(wls)
% This function calculates the diffraction angle after passing through the
% dwp module in radians.
% 
% We use the Sellmeier equation to calculate the refractive performance of
% the glass wrt to different wavelengths.
%
% Constants used are based on the DWP module specfications. Prism angle =
% 30 deg. Glass1 = N-SF11. Glass2 = N-LAF21.
%
% Input:
%   wls: wavelength in um
%
% Output:
%   pos: diffraction angle in radians
%
%
% Created by Weihong Yeo, Northwestern University, 2021-12-29.
% Last modified by Weihong Yeo, Northwestern University, 2023-05-18.
%


prism_angle = 30;

eta0 = 1.0003;
eta1 = sellmeier(wls, 'N-SF11');
eta2 = sellmeier(wls, 'N-LAF21');

% interface between first and second prism
toprism2_angle = asind(eta1./eta2 .* sind(prism_angle));

% interface between second prism and air
pos = asin(eta2./eta0 .* sind(prism_angle-toprism2_angle));

end


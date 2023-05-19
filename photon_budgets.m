function [B0,B1] = photon_budgets(photons, spectrometer, debug)
% This function calculates the background photon budget for a given splitting
% ratio.
%
% Inputs:
%   photons: a struct containing the number of photons in each channel
%   spectrometer: a struct containing the spectrometer information
%   debug: a boolean to plot the efficiency functions
%
% Outputs:
%   B0: the number of background photons in channel 0
%   B1: the number of background photons in channel 1
%
%
% Created by Weihong Yeo, Northwestern University, 2021-12-29.
% Last modified by Weihong Yeo, Northwestern University, 2023-05-18.
%

if nargin == 2
    debug = false;
end
if photons.B == 0
    B0 = 0;
    B1 = 0;
    
elseif photons.B_equal
    B0 = photons.B / 2;
    B1 = photons.B / 2;
    
else
%% 
% efficiency functions should have been loaded.
    x = spectrometer.wl_range(1):spectrometer.wl_range(2);
    y0 = spectrometer.order0_fit(x);
    y1 = 1 - spectrometer.order1_fit(x);
    
    if debug == true
        plot(x,y0,x,y1);
    end
    
    B0 = photons.B * trapz(x,y0) / (spectrometer.wl_range(2) - spectrometer.wl_range(1));
    B1 = photons.B * trapz(x,y1) / (spectrometer.wl_range(2) - spectrometer.wl_range(1));
    
end
end
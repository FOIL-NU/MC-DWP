function [out, uncaptured] = psfsample(x,y,s,px_sz)
% This function implements the sampling of the PSF onto the imaging sensor
% for a given set of photon positions. The PSF is assumed to be a
% gaussian, and the sampling is done by discretizing the PSF into pixels
% and then counting the number of photons that land on each pixel.
%
% Note:
% we cannot have an ideal psf to sample because the splitting efficiency is
% different at different wavelengths (for the grating case), therefore we
% need do a different approach whereby we treat each photon individually.
%
% Inputs:
%   x: vector of x positions of photons landing on the imaging sensor
%   y: vector of y positions of photons landing on the imaging sensor
%   s: size of the image sensor crop region
%   px_sz: pixel size
%
% Outputs:
%   out: sampled PSF on the camera sensor
%   uncaptured: number of photons that did not land on the imaging sensor
%
%
% Created by Weihong Yeo, Northwestern University, 2021-12-29.
% Last modified by Weihong Yeo, Northwestern University, 2023-05-18.
%

% list full size
sx = s(1); sy = s(2);

% psf_edges
psf_edgesx = [-Inf, px_sz*((0:sx)-sx/2), Inf];
psf_edgesy = [-Inf, px_sz*((0:sy)-sy/2), Inf];

% get locations of each photon on the pixel array
% mod changes out of range pixels to zero.
loc_x = mod(discretize(x,psf_edgesx)-1,sx+1);
loc_y = mod(discretize(y,psf_edgesy)-1,sy+1);

% if the pixel is out of range, set loc_x = 0 and loc_y = 1
loc_x(loc_y == 0) = 1;
loc_y(loc_x == 0) = 0;

% get the count of each photon on each pixel, then reshape
% if loc_x = 0 and loc_y = 1, histcounts accumulates to 0.
locs = histcounts(loc_x+(loc_y-1)*sx,'BinMethod','integers','BinLimits',[0,prod(s)]);

% reshape everything except for the out of range pixels.
out = reshape(locs(2:end),s)';

% optional output to return the number of wasted (uncaptured) photons.
uncaptured = locs(1);
if uncaptured/numel(x) > 0.01
    warning('%.1f%% of photons are uncaptured. Consider changing crop boundary!', ...
        uncaptured/numel(x)*100);
end
function outdata = ssmlm_sim( ...
        photons, sim_settings, specmtr, spectral, camera, microscope, n_iter, rng_init)
% This function simulates the sSMLM experiment. It takes in the parameters
% of the experiment, and outputs the calculated parameters of the experiment.
%
% Inputs:
%   photons: struct containing the parameters of the photons
%       .N: number of photons
%       .mode: 'lognormal' or 'constant'
%       .N_mu: mean of lognormal distribution
%       .N_sig: standard deviation of lognormal distribution
%       .bg: background photons
%       .bg_per_px: background photons per pixel
%       .bg_per_px_0: background photons per pixel for order 0
%       .bg_per_px_1: background photons per pixel for order 1
%   sim_settings: struct containing the parameters of the simulation
%       .psf_size: size of the psf in pixels
%   specmtr: struct containing the parameters of the spectrometer
%       .type: 'dwp' or 'grating'
%       .specdispersion: spectral dispersion of the spectrometer
%   spectral: struct containing the parameters of the fluorophore
%       .type: 'nanodrop', 'nilered', 'afxxx', 'cyxxx', 'cfxxx', 'dyxxx',
%              'tritc', 'gaussian', 'hetero', 'hetereo_gaussian'
%       .lambda: peak wavelength of the fluorophore
%       .sigma: width of the fluorophore
%       .cov: covariance matrix of the fluorophore if gaussian model is used
%       .det_min: minimum wavelength of the spectrometer
%       .det_max: maximum wavelength of the spectrometer
%   camera: struct containing the parameters of the camera
%       .actualpxsize: actual pixel size of the camera
%       .magnification: magnification of the microscope
%       .gain: gain of the camera
%       .adu: adu of the camera
%   microscope: struct containing the parameters of the microscope
%       .NA: numerical aperture of the microscope
%   n_iter: number of iterations
%   rng_init: random number generator seed
%
% Outputs:
%   outdata: struct containing the parameters of the experiment
%       .gof0: goodness of fit for order 0
%       .gof1: goodness of fit for order 1
%       .coeff0: coefficients for the gaussian fit of order 0
%       .coeff1: coefficients for the gaussian fit of order 1
%       .photons_n: number of photons
%       .order0_raw: raw image of order 0
%       .order1_raw: raw image of order 1
%       .order1_wls: corresponding wavelengths of order 1
%       .order0_peak_estd: estimated peak location of order 0
%       .order1_peak_estd: estimated peak location of order 1
%       .order0_width_estd: estimated width of order 0
%       .order1_width_estd: estimated width of order 1
%       .order0_peak_true: true peak location of order 0
%       .order1_peak_true: true peak location of order 1
%       .order0_width_true: true width of order 0
%       .order1_width_true: true width of order 1
%       .order0_loc_prec: localization precision of order 0
%       .order1_centroid_estd: estimated centroid of order 1 using spectral 
%           centroid
%
%
% Created by Weihong Yeo, Northwestern University, 2021-12-29.
% Last modified by Weihong Yeo, Northwestern University, 2023-05-18.
%
% #########################################################################
% Changelog
% #########################################################################
% 
%  2022-01-17
%   - added gaussian fluorophore spectral heterogeneity input (change in peak)
%  
%  2022-02-08
%   - changed spectral heterogeneity input to covariance matrix for peaks and widths
%   - added varying n_photon counts based on lognormal distribution
%   - cleaned output parameters: removed useless output parameters, and renamed parameters
%  
%  2022-03-22
%   - migrated camera.magnification to microscope.magnification
%   - added variables for gain and adu for camera
%  
%  2022-04-27
%   - per-photon psf used instead of an average value
%  
%  2022-09-09
%  - corrected spectral calculation. scales the intensity based on spectral dispersion now.
%

% check inputs
if nargin == 8
    rng(rng_init);
end

if length(sim_settings.psf_size) == 1
    sim_settings.psf_size = repmat(sim_settings.psf_size,1,2);
elseif length(sim_settings.psf_size) > 2
    error('psf_size in sim_settings is incorrect!');
end

% check if savefile exists
if strcmp(photons.mode, 'lognormal')
    photons.N = round(exp(photons.N_mu + photons.N_sig ^2 / 2) * 2);
end

[file_exists, file_path] = checksavefile(photons, sim_settings, specmtr, spectral, camera, microscope, n_iter, rng_init);

if file_exists
    load(file_path);
else
% evaluate parameters
mean_psf_lambda = spectral.lambda;
mean_psf_width = 0.61 * mean_psf_lambda / microscope.NA;
mean_psf_sigma = mean_psf_width / 2.355;
camera_pxsize = camera.actualpxsize / microscope.magnification;
% load fluorophore data
% if spectral hetereogeneity is selected, then we resample the spectral distribution 
% every iteration
spectral_dist = cell(n_iter, 1);
order1_peak_true = ones(n_iter,1) * spectral.lambda;
order1_width_true = ones(n_iter,1) * spectral.sigma;
if contains(spectral.type, 'hetereo')
    disp('spectral hetereogenity');
    for iter = 1:n_iter
        spectral = load_spectral(spectral);
        spectral_dist{iter} = spectral.dist;
        order1_peak_true(iter) = spectral.dist.mu;
        order1_width_true(iter) = spectral.dist.sigma;
    end
else
    spectral = load_spectral(spectral);
    for iter = 1:n_iter
        spectral_dist{iter} = spectral.dist;
    end
end

% load grating/dwp module efficiency data
specmtr = load_efficiency(specmtr);
%% 
% generate centroid of the fluorescence
pixel_centroids_true = random('Uniform',-camera_pxsize/2,camera_pxsize/2,n_iter,2);
% dwp module
% for some reason it wouldnt run if this was not loaded, probably because it 
% is needed when evaluating dwp fit
white_light = linspace(spectral.det_min,spectral.det_max)/1000;
out_angle = dwp_module(white_light);
%% 
% calculate the distance away from the camera required for the target spectral 
% dispersion
dwp_pos = (16e-6 * (white_light(end)-white_light(1)) * 1000) / ...
    (specmtr.specdispersion * .01 * (tan(out_angle(end) - out_angle(1))));
dwp_spread = .01 * dwp_pos / 16e-6 * tan(out_angle);
%% 
% adjust back the dwp_spread so that all the photons are centered
dwp_spread = dwp_spread - (dwp_spread(1) + dwp_spread(end)) / 2;
%% 
% correlate white light and dwp_spread
dwp_fit = fit(white_light(:)*1000, dwp_spread(:), 'smoothingspline');
dwp_invfit = fit(dwp_spread(:), white_light(:)*1000, 'smoothingspline');
% change variables stored in struct into normal variables
if strcmp(photons.mode, 'lognormal')
    photons_n = round(random('lognormal', photons.N_mu, photons.N_sig, n_iter, 1) * 2);
else
    photons_n = photons.N * ones(n_iter,1);
end
[bg_per_px_0, bg_per_px_1] = photon_budgets(photons,specmtr);

% change variables stored in struct into normal variables so that we can
% use parfor.
sim_settings_psf_size = sim_settings.psf_size;
specmtr_type = specmtr.type;
specmtr_bs_split = specmtr.bs_split;
specmtr_sd = specmtr.specdispersion;
specmtr_order0_fit = specmtr.order0_fit;
specmtr_order1_fit = specmtr.order1_fit;
spectral_lambda = spectral.lambda;
spectral_sigma = spectral.sigma;
spectral_det_min = spectral.det_min;
spectral_det_max = spectral.det_max;

microscope_NA = microscope.NA;
camera_readout_sigma = camera.readout_sigma;
camera_gain = camera.gain;
camera_adu = camera.adu;

% camera spectral cutoff
% calculate spectral cutoff so that we can correlate the image with a wavelength
if contains(specmtr_type, 'grating')
    % all in units of wavelength (nm)
    camera_spectralpx_edges = spectral_det_min:specmtr_sd:spectral_det_max;
    camera_spectralpx_center = (camera_spectralpx_edges(2:end) + camera_spectralpx_edges(1:end-1)) / 2;
    camera_spectralpx_dispersion = diff(camera_spectralpx_edges);
elseif contains(specmtr_type, 'dwp')
    % indices units of pixels
    camera_spectralpx_indices = floor(min(dwp_spread)):(ceil(max(dwp_spread))-1)+0.5;
    % center in units of wavelength (nm)
    camera_spectralpx_center = dwp_invfit(camera_spectralpx_indices)';
    camera_spectralpx_dispersion = diff(dwp_invfit(floor(min(dwp_spread)):(ceil(max(dwp_spread)))))';
end
% spectral dispersion size
% calculate size of simulated localization
psf1_wl_width = length(camera_spectralpx_center) - 1;
psf0_size = sim_settings_psf_size*2+1;
psf1_size = sim_settings_psf_size*2+1 + [psf1_wl_width,0];
sim_settings_psf_size_1 = sim_settings_psf_size(1);
sim_settings_psf_size_2 = sim_settings_psf_size(2);
%% 
% arrays to store data
spectral_centroids = zeros(n_iter,1);
gof0 = nan(n_iter,1);
gof1 = nan(n_iter,1);
coeff0 = nan(n_iter,5);
coeff1 = nan(n_iter,6);
roi0 = nan([n_iter,psf0_size]);
roi1 = nan([n_iter,psf1_size]);
wls = nan(n_iter,psf1_wl_width+1);

% start simulation
parfor iter = 1:n_iter
% spectral distribution
% generate wavelength of each of the _N_ photons based on the specified distribution. 
% remove photons that are out of the concerned range
    photon_wl = random(spectral_dist{iter},photons_n(iter),1);
%% 
% photon-level psf
    psf_width = 0.61 * photon_wl / microscope_NA;
    psf_sigma = psf_width / 2.355;
%% 
% generates the location of each of the _N_ photons sampled based on a gaussian 
% psf
    photon_psfshift = pixel_centroids_true(iter,:)+psf_sigma.*randn(photons_n(iter),2);
%% 
% deflects the photon based on the type of spectrometer specified
    if contains(specmtr_type, 'grating')
        spectral_det_avg = (spectral_det_max + spectral_det_min) / 2;
        photon_wldeflect = (photon_wl-spectral_det_avg)/specmtr_sd*camera_pxsize;
    elseif contains(specmtr_type, 'dwp')
        photon_wldeflect = feval(dwp_fit,photon_wl)*camera_pxsize;
    end
%% 
% calculate photon split order
    photon_order0_cutoff = specmtr_order0_fit(photon_wl);
    photon_order1_cutoff = specmtr_order1_fit(photon_wl);
    photon_order_rv = rand(photons_n(iter),1);
%% 
% order_rv < order0_cutoff, set to order0 (0), order1_cutoff < order_rv, set 
% to order1 (1), else lost (-1).
    photon_order = -ones(photons_n(iter), 1); % initialize with -1, all lost
    photon_order(photon_order_rv < photon_order0_cutoff) = 0;
    photon_order(photon_order1_cutoff < photon_order_rv) = 1;
    photon_order((photon_wl<spectral_det_min) & (photon_wl>spectral_det_max)) = -1;
% generating the image
% signal noises
% shot noise sample of psf0 and psf1
    psf0_data = psfsample( ...
        photon_psfshift(photon_order == 0,1), ...
        photon_psfshift(photon_order == 0,2), ...
        psf0_size,camera_pxsize);
    psf1_data = psfsample( ...
        photon_psfshift(photon_order == 1,1) + photon_wldeflect(photon_order == 1), ...
        photon_psfshift(photon_order == 1,2), ...
        psf1_size,camera_pxsize);
%% 
% EMCCD amplification
    psf0_data = random('Gamma',psf0_data,camera_gain)/camera_adu;
    psf1_data = random('Gamma',psf1_data,camera_gain)/camera_adu;
%% 
% readout noise
    psf0_data = random('Normal',0,camera_readout_sigma,size(psf0_data)) + psf0_data;
    psf1_data = random('Normal',0,camera_readout_sigma,size(psf1_data)) + psf1_data;
% background noise
    if (bg_per_px_0 > 0) && (bg_per_px_1 > 0)
%% 
% Poisson noise for background, and adding the same EMCCD amplification
        psf0_bg = random('Poisson',bg_per_px_0,size(psf0_data));
        psf0_bg = random('Gamma',psf0_bg,camera_gain)/camera_adu;
        
        psf1_bg = random('Poisson',bg_per_px_1,size(psf1_data));
        psf1_bg = random('Gamma',psf1_bg,camera_gain)/camera_adu;
%% 
% remove expected background noise
        psf0_bg = psf0_bg - bg_per_px_0*ones(size(psf0_data))*camera_gain/camera_adu;
        psf1_bg = psf1_bg - bg_per_px_1*ones(size(psf1_data))*camera_gain/camera_adu;
%% 
% add the background to the raw camera sensor data
        psf0_data = psf0_data + psf0_bg;
        psf1_data = psf1_data + psf1_bg;
    end
%% 
% round off because all values are integral
    psf0_data = round(psf0_data);
    psf1_data = round(psf1_data);
    roi0(iter,:,:) = psf0_data';
    roi1(iter,:,:) = psf1_data';
%% 
% compute the ysum of the sensor data
%     psf0_ysum = sum(psf0_data,1);
    psf1_ysum = sum(psf1_data,1);
    
    % crop the relevant regions out, otherwise it will be too long
    psf1_ysum = psf1_ysum(sim_settings_psf_size_1+1:end-sim_settings_psf_size_1);
    psf1_crop = psf1_data(:,sim_settings_psf_size_1+1:end-sim_settings_psf_size_1);
% estimate pixel centroids
% perform 2d gaussian fit to estimate the pixel centroids
    if ~strcmpi(specmtr_bs_split,'none')
        x0 = -sim_settings_psf_size_1:sim_settings_psf_size_1;
        y0 = -sim_settings_psf_size_2:sim_settings_psf_size_2;
        
        f0_upper = [1, 1, 1, sim_settings_psf_size_1/2, 5];
        f0_lower = [0, -1, -1, 0, 0];
        f0_initguess = [1/(2*pi*(mean_psf_sigma/camera_pxsize)^2), ...
            0, 0, mean_psf_sigma/camera_pxsize, 0];
        f0_X = [repelem(x0(:),length(y0)), repmat(y0(:),length(x0),1)];
        f0_Y = psf0_data(:);
        f0_Y(f0_Y < 0) = 0;
        f0_scale = sum(f0_Y);
        f0_Y = f0_Y ./ f0_scale;
        
        fn0 = @(params)(mle_fn(f0_X,f0_Y,params));
        vals = fminsearch(fn0,f0_initguess,optimset('display','none'));
        coeff0(iter,:) = [f0_scale, 1, 1, 1, 1] .* vals;
        pixel_centroids_estd = vals(2);

    else
        pixel_centroids_estd = 0;
        
    end
% estimate spectral centroids
    x1 = camera_spectralpx_center-pixel_centroids_estd/camera_pxsize*specmtr_sd;
    y1 = -sim_settings_psf_size_2:sim_settings_psf_size_2;
    
    wls(iter,:) = x1;
    
    spectral_centroids(iter) = sum(psf1_ysum .* x1)/sum(psf1_ysum);
    % we perform least squares fit to estimate the peak instead of mle because it is
    % faster and more robust
    % we use the spectral centroid as the initial guess for the spectral
    f1_upper = [photons_n(iter), spectral_det_max, 2, ...
        (spectral_det_max - spectral_det_min)/4, sim_settings_psf_size_2/2, ...
        5];
    f1_lower = [0, spectral_det_min, -2, 0, 0, -5];
    f1_initguess = [max(psf1_data,[],'all'), spectral_lambda, pixel_centroids_true(iter,2), ...
        spectral_sigma, mean_psf_sigma/camera_pxsize, 0];
    f1_X = [repelem(x1(:),length(y1)), repmat(y1(:),length(x1),1)];
    f1_Y = psf1_crop(:);
    f1_Y = f1_Y./repelem(camera_spectralpx_dispersion(:),length(y1));
    fn1 = fittype('gaussian2dellipse(X,Y,a,b1,b2,c1,c2,d)', ...
        'independent', {'X','Y'}, ...
        'coefficients', {'a','b1','b2','c1','c2','d'});
    try
        [f1,gof1_] = fit(f1_X,f1_Y,fn1, ...
            'StartPoint',f1_initguess, ...
            'Lower', f1_lower, ...
            'Upper', f1_upper, ...
            'MaxIter', 1000, ...
            'TolFun', 1e-15, ...
            'TolX', 1e-10);
        
        coeff1(iter,:) = coeffvalues(f1);
        gof1(iter) = gof1_.rsquare;
        
    catch
        % do nothing, NaN assumed in the array creation.
        
    end
end % parfor
% accumulate the output data into the structure |outdata|
outdata.gof0 = gof0; % used for old code compatibility
outdata.gof1 = gof1;
outdata.coeff0 = coeff0;
outdata.coeff1 = coeff1;
outdata.photons_n = photons_n;
outdata.order0_raw = roi0;
outdata.order1_raw = roi1;
outdata.order1_wls = wls;
%% 
% for reference/compatibility, but note that they are mostly renamed
outdata.order0_peak_estd = coeff0(:,2:3);
outdata.order0_width_estd = coeff0(:,4:5);
outdata.order0_peak_true = pixel_centroids_true;
outdata.order0_width_true = mean_psf_sigma;
outdata.order0_loc_prec = mean(std((coeff0(:,2:3).*camera_pxsize - pixel_centroids_true),1));
%% 
% where localization precision is the mean of x and y
outdata.order1_peak_estd = coeff1(:,2:3);
outdata.order1_width_estd = coeff1(:,4:5);
outdata.order1_peak_true = order1_peak_true;
outdata.order1_width_true = order1_width_true;
outdata.order1_centroid_estd = spectral_centroids;
%% 
% note that the means width is a convolution of the psf of the zeroth order 
% and the first order spectral width. and orig, estd, true should have similar 
% result for the specwidth calculation, because orig, estd, true only refers to 
% the sse correction used.

end
save(file_path);
end
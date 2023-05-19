function spectral = load_spectral(spectral, debug)
% This function loads the spectral information from the fpbase database,
% and returns a structure with the spectral information. We also allow for
% a gaussian distribution to be used, or custom nanodrop data to be used.
%
% Inputs:
%   spectral: struct containing the parameters of the fluorophore
%       .type: 'nanodrop', 'nilered', 'afxxx', 'cyxxx', 'cfxxx', 'dyxxx',
%              'tritc', 'gaussian', 'hetero', 'hetereo_gaussian'
%       .lambda: peak wavelength of the fluorophore
%       .sigma: width of the fluorophore
%       .cov: covariance matrix for the peak and width of the fluorophore
%   debug: boolean to plot the spectral distribution
%
% Outputs:
%   spectral: appended struct containing the parameters of the fluorophore
%       .data: struct containing the data from the fpbase database
%       .dist: distribution of the fluorophore
%
% Created by Weihong Yeo, Northwestern University, 2021-12-29.
% Last modified by Weihong Yeo, Northwestern University, 2023-05-18.
%
% #########################################################################
% Changelog
% #########################################################################
% 
%  2021-12-29
%   - initial creation
%
%  2022-01-17 
%   - add gaussian fluorophore spectral heterogeneity input (change in peak)
% 
%  2022-02-08
%   - change spectral heterogeneity input to covariance matrix for peaks and widths
%
%  2022-05-24 
%   - added nanodrop processing code
%

% check inputs
if nargin == 1
    debug = false;
end
% nanodrop data
if contains(spectral.type,'nanodrop')
    root_folder = 'nanodrop';
    if contains(spectral.type,'ns')
        spectral.data = load_fluor(sprintf('%s/ns%s.csv', root_folder, spectral.type(end-2:end)));
        spectral.dist = fitdist(spectral.data.wl, 'Kernel', 'Width', 2, 'frequency', round(spectral.data.em*1000));
    end
else
    root_folder = 'fpbase';
% fpbase data
% nile red
    if contains(spectral.type, 'nilered')
        if contains(spectral.type, 'p')
            spectral.data = load_fluor(sprintf('%s/Nile Red (phospholipid).csv', root_folder));
            
        elseif contains(spectral.type, 't')
            spectral.data = load_fluor(sprintf('%s/Nile Red (triglyceride).csv', root_folder));
        
        else
            error('not implemented, expected either (t)riglyceride or (p)hospholipid');
            
        end
        
        spectral.dist = fitdist(spectral.data.wl, 'Kernel', 'Width', 1, 'frequency', round(spectral.data.em*1000));
%% 
% alexa fluor series
    elseif contains(spectral.type, 'af')
        spectral.data = load_fluor(sprintf('%s/Alexa Fluor %s.csv', root_folder, spectral.type(end-2:end)));
        spectral.dist = fitdist(spectral.data.wl, 'Kernel', 'Width', 2, 'frequency', round(spectral.data.em*1000));
%% 
% cy series
    elseif contains(spectral.type, 'cy')
        spectral.data = load_fluor(sprintf('%s/Cy%s.csv', root_folder, spectral.type(3:end)));
        spectral.dist = fitdist(spectral.data.wl, 'Kernel', 'Width', 1, 'frequency', round(spectral.data.em*1000));
%% 
% cf series
    elseif contains(spectral.type, 'cf')
        spectral.data = load_fluor(sprintf('%s/%s.csv', root_folder, upper(spectral.type)));
        spectral.dist = fitdist(spectral.data.wl, 'Kernel', 'Width', 1, 'frequency', round(spectral.data.em*1000));   
%% 
% dy series
    elseif contains(spectral.type, 'dy')
        spectral.data = load_fluor(sprintf('%s/DY-%s.csv', root_folder, spectral.type(end-2:end)));
        spectral.dist = fitdist(spectral.data.wl, 'Kernel', 'Width', 1, 'frequency', round(spectral.data.em*1000));
%% 
% tritc
    elseif contains(spectral.type, 'tritc')
        spectral.data = load_fluor(sprintf('%s/tritc.csv', root_folder));
        spectral.dist = fitdist(spectral.data.wl, 'Kernel', 'Width', 1, 'frequency', round(spectral.data.em*1000));    
%% 
% gaussian
    elseif contains(spectral.type, 'gaussian')
        if contains(spectral.type, 'hetereo')
            spectral_params = zeros(1,2);
            while any(spectral_params <= 0)
                spectral_params = mvnrnd([spectral.lambda, spectral.sigma], spectral.cov);
            end
            spectral.dist = makedist('Normal', spectral_params(1), spectral_params(2));
        else
            spectral.dist = makedist('Normal', spectral.lambda, spectral.sigma);
        end
%% 
% catch the rest
    elseif isfile(sprintf('%s/%s.csv', root_folder, spectral.type))
        spectral.data = load_fluor(sprintf('%s/%s.csv', root_folder, spectral.type));
        spectral.dist = fitdist(spectral.data.wl, 'Kernel', 'Width', 1, 'frequency', round(spectral.data.em*1000));
    else
        error('not implemented, expected either afxxx, cyx, nilered or gaussian, or file not found.');
    end
%% 
% plot graphs for debugging
if debug
    wl = spectral.det_min:spectral.det_max;
    em = pdf(spectral.dist, wl);
    
    plot(wl, em);
    xlabel('wavelength (nm)');
    ylabel('emission (arb)');
end
end
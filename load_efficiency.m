function spectrometer = load_efficiency(spectrometer, debug)
% This function loads the efficiency data for the spectrometer. The
% efficiency data is stored in the spectrometer structure as a fit object.
%
% Inputs:
%   spectrometer: structure containing the spectrometer information.
%   debug: boolean, whether to plot the efficiency data.
%
% Outputs:
%   spectrometer: structure containing the spectrometer information, with
%   the efficiency data stored as a fit object.
%
%
% Created by Weihong Yeo, Northwestern University, 2022-03-23.
% Last modified by Weihong Yeo, Northwestern University, 2022-03-23.
%

% check inputs
if nargin == 1
    debug = false;
end

% load custom ideal efficiency curve
if contains(spectrometer.type, 'ideal')
    wl = spectrometer.wl_range;
    eff0 = spectrometer.ideal_params.eff0;
    eff1 = 1 - spectrometer.ideal_params.eff1;
    fit_type = 'poly1';
    disp('ideal efficiency');
    
    % correct the inputs of the fits.
    if length(wl) == 1
        error('wl should be a 1x2 vector for the range of wavelengths');
        
    else
        if length(eff0) == 1
            eff0 = repmat(eff0, size(wl));
        end
        
        if length(eff1) == 1
            eff1 = repmat(eff1, size(wl));
        end
    end
%% 
% efficiency curve for staranalyzer100
elseif contains(spectrometer.type, 'grating')
    load('grating_data.mat', 'grating_data');
    disp('real grating efficiency loaded')
    
    wl = grating_data(:,1);
    eff0 = grating_data(:,2);
    eff1 = 1 - grating_data(:,3);
    fit_type = 'poly1';
%% 
% efficiency curve for dwp module
elseif strcmp(spectrometer.type, 'dwp')
    available_splits = {'1090', '3070', '5050', '7030', '9010'};
    split = spectrometer.bs_split;
    wl_range = spectrometer.wl_range;
    
%% 
% load dwp and ra_prism efficiency data
    eff_data = readtable('eff_data.csv','PreserveVariableNames',true);
    wl = eff_data{:, 'wavelength (nm)'};
%% 
% for dwp, smooth data with sgolay filter, then perform interpolation within 
% the range. for values beyond 690nm, the last value is used to extrapolate.
    dwp_eff = smoothdata(eff_data{:, 'dwp'},'sgolay');
    dwp_eff = interp1(wl,dwp_eff,wl_range(1):wl_range(2),'pchip',dwp_eff(end));
%% 
% for right angle prism, take the mean of the efficiency values and assume it 
% is constant.
    rap_eff = mean(eff_data{:, 'ra_prism'});
%% 
% load beamsplitter efficiency data
    if any(strcmpi(split,available_splits))
        rawdata = readtable('bs_data.csv', 'PreserveVariableNames', true);
        wl = rawdata{:, 'wavelength (nm)'};
        relevant_wl = (wl>=wl_range(1)) & (wl<=wl_range(2));
        p_transmit = rawdata{relevant_wl, strcat(split, ' p transmit')};
        p_reflect = rawdata{relevant_wl, strcat(split, ' p reflect')};
        s_transmit = rawdata{relevant_wl, strcat(split, ' s transmit')};
        s_reflect = rawdata{relevant_wl, strcat(split, ' s reflect')};
        
        % first order
        mean_transmit = 1 - mean([p_transmit, s_transmit], 2)/100 .* dwp_eff(:);
        % zeroth order
        mean_reflect = mean([p_reflect, s_reflect], 2)/100 .* rap_eff(:);
        
        wl = wl_range(1):wl_range(2);
        eff0 = mean_reflect;
        eff1 = mean_transmit;
        fit_type = 'smoothingspline';
        
    elseif strcmpi(split,'none')
        wl = wl_range(1):wl_range(2);
        eff0 = zeros(1,length(dwp_eff));
        eff1 = 1 - dwp_eff;
        fit_type = 'smoothingspline';
        
    else
        error('invalid split! please input "%s", for example.', available_splits{1});
        
    end
else
    error('not implemented.');
end
spectrometer.order0_fit = fit(wl(:), eff0(:), fit_type);
spectrometer.order1_fit = fit(wl(:), eff1(:), fit_type);
%% 
% plot graphs for debugging
if debug
    if contains(specmtr.type, 'grating')
        figure(1);
        clf;
        hold on;
        scatter(specmtr.data_wl, specmtr.data_eff0);
        scatter(specmtr.data_wl, specmtr.data_eff1);
        plot(specmtr.data_wl, polyval(specmtr.order0_fit, specmtr.data_effwl));
        plot(specmtr.data_wl, polyval(specmtr.order1_fit, specmtr.data_effwl));
    
    end
end
end
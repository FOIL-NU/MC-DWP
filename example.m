%% sample script
%% experimental parameters

n_iterations = 100;
photons.N = 1000;  % photon budget for both channels
photons.B = 10;    % photons per pixel for both channels
photons.B_equal = false;
photons.mode = 'constant';

% place varied parameters below...
parameter_in1 = [680];

%% default parameters
% spectral characteristics
spectral.type = 'gaussian'; % 'gaussian', or just the fluorophore name
spectral.lambda = 670; % nm
spectral.sigma  = 20;  % nm
spectral.cov    = zeros(2);

spectral.det_min = 450; % nm
spectral.det_max = 850; % nm
% camera characteristics

camera.actualpxsize = 16000; % nm
camera.readout_sigma = 2;
camera.gain = 100;
camera.adu = 14.3;
% spectrometer characteristics

spectrometer.type = 'dwp';

spectrometer.bs_split = '5050';
spectrometer.wl_range = [spectral.det_min,spectral.det_max];

spectrometer.specdispersion = 4.3; % nm/px

spectrometer.ideal_params.eff0 = 0.265;
spectrometer.ideal_params.eff1 = 0.550;
% microscope characteristics

microscope.NA = 1.49;
microscope.magnification = 100;
% simulation settings

sim_settings.psf_size = [4,4];
sim_settings.order0_fitmethod = 'average';
sim_settings.order1_fitmethod = 'average';
%% run simulation
% prepare output parameters

parameter_out = cell(1,length(parameter_in1));
%% 
% overwrite the relevant parameters in the loop

for i1 = 1:length(parameter_in1)
    spectral.lambda = parameter_in1(i1);
    
    spectrometer.type = 'dwp';
    spectrometer.bs_split = 'none';
    
    outdata = ssmlm_sim(photons, sim_settings, ...
        spectrometer, spectral, camera, microscope, n_iterations, 8);
    parameter_out{i1} = outdata;
end
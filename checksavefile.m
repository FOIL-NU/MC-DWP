function [file_exists, file_path] = checksavefile( ...
    photons, sim_settings, specmtr, spectral, camera, microscope, n_iterations, rng_init)
% This function checks if a savefile exists for the current simulation
% settings. If it does, it returns the path to the file. If it does not, it
% returns the path to the next file to be created.
% Stores savefiles in a structure: 
%   data/iterxxx/BSxxxx/Bxx/Nxx/lambdaxxx/widthxxx/001.mat
%
% Inputs:
%   same as ssmlm_sim.m
%
% Outputs:
%   file_exists: boolean, true if file exists, false if not
%   file_path: path to the file, whether it exists or not
%
%
% Created by Weihong Yeo, Northwestern University, 2021-12-29.
% Last modified by Weihong Yeo, Northwestern University, 2023-05-18.
%

currdir = 'data';
chkmkdir(currdir);

currdir = [currdir, '/iter', sprintf('%d', n_iterations)];
chkmkdir(currdir);

if contains(specmtr.type,'grating')
    currdir = [currdir, '/grating'];
else
    currdir = [currdir, '/BS', specmtr.bs_split];
end
chkmkdir(currdir);

if photons.B_equal
    currdir = [currdir, '/B', sprintf('%d',photons.B), 'equal'];
else
    currdir = [currdir, '/B', sprintf('%d',photons.B)];
end
chkmkdir(currdir);

currdir = [currdir, '/N', sprintf('%d',photons.N)];
chkmkdir(currdir);

currdir = [currdir, '/lambda', sprintf('%d',round(spectral.lambda))];
chkmkdir(currdir);

currdir = [currdir, '/width', sprintf('%d',round(spectral.sigma))];
chkmkdir(currdir);

file_index = 1;
file_exists = false;
file_path = [currdir, '/', sprintf('%03d.mat', file_index)];

while ~file_exists && isfile(file_path)
    existing_file = load(file_path);
    
    equal_sim_settings = isequal(existing_file.sim_settings, sim_settings);
    
    % check specmtr settings
    equal_specmtr = true;
    
    if ~isequal(existing_file.specmtr.type, specmtr.type)
        equal_specmtr = false;
    elseif ~isequal(existing_file.specmtr.bs_split, specmtr.bs_split)
        equal_specmtr = false;
    elseif ~isequal(existing_file.specmtr.wl_range, specmtr.wl_range)
        equal_specmtr = false;
    elseif ~isequal(existing_file.specmtr.specdispersion, specmtr.specdispersion)
        equal_specmtr = false;
    elseif ~isequal(existing_file.specmtr.ideal_params, specmtr.ideal_params)
        equal_specmtr = false;
    end
    
    % check spectral settings
    equal_spectral = true;
    
    if ~isequal(existing_file.spectral.cov, spectral.cov)
        equal_spectral = false;
    end
    
    if ~isequal(existing_file.spectral.type, spectral.type)
        equal_spectral = false;
        disp('unequal spectral type');
    elseif ~isequal(existing_file.spectral.lambda, spectral.lambda)
        equal_spectral = false;
        disp('unequal spectral lambda');
    elseif ~isequal(existing_file.spectral.sigma, spectral.sigma)
        equal_spectral = false;
        disp('unequal spectral sigma');
    elseif ~isequal(existing_file.spectral.det_min, spectral.det_min)
        equal_spectral = false;
        disp('unequal spectral det_min');
    elseif ~isequal(existing_file.spectral.det_max, spectral.det_max)
        equal_spectral = false;
        disp('unequal spectral det_max');
    end

    
    equal_camera = isequal(existing_file.camera, camera);
    if ~equal_camera
        disp('unequal camera');
    end
    equal_microscope = isequal(existing_file.microscope, microscope);
    if ~equal_microscope
        disp('unequal microscope');
    end
    equal_rng_init = isequal(existing_file.rng_init, rng_init);
    if ~equal_rng_init
        disp('unequal rng');
    end
    equal_checks = [equal_sim_settings, equal_specmtr, equal_spectral, ...
        equal_camera, equal_microscope, equal_rng_init];
    
    if all(equal_checks)  % if all the data are equal, 
        file_exists = true;
        file_path = [currdir, '/', sprintf('%03d.mat', file_index)];
        fprintf('file exists, loaded from %s.\n', file_path);
    else
        file_index = file_index + 1;
        file_path = [currdir, '/', sprintf('%03d.mat', file_index)];
    end
end

if ~file_exists
    fprintf('file not found, creating %s.\n', file_path);
end

end

function chkmkdir(dir)
    if ~isfolder(dir)
        mkdir(dir)
    end
end

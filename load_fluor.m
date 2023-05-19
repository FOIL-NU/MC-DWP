function fluor = load_fluor(filename)
% Helper function to load fluorophore data from a file
%
% Inputs:
%   filename: path to file containing fluorophore data
%
% Outputs:
%   fluor: struct containing fluorophore data
%      .wl: wavelengths (nm)
%      .em: emission intensity in arbitrary units
%      .ex: excitation intensity in arbitrary units
%
%
% Created by Weihong Yeo, Northwestern University, 2022-03-23.
% Last modified by Weihong Yeo, Northwestern University, 2022-03-23.
%

%% read fluorophore data
out = readtable(filename, 'PreserveVariableNames', true);

isnonzero = (out{:,'Emission'} >= 0.01);

fluor.wl = out{isnonzero,'Wavelength'};
fluor.em = out{isnonzero,'Emission'};
fluor.ex = out{isnonzero,'Excitation'};

end
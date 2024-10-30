%%% Set filepath for the functions
clear; close all; fig = 0;

disp('!Attention!');
disp("Run the main code when inside the 'scripts' folder");
disp('so that the filepaths are setup correctly.');
disp(char(9) + "Thank you."); disp(' ');

mainPath = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(mainPath, '..', 'src')));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load files - mri t1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading files...');

load_files              % Loads MRI scan    

disp(' ');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the final 3D model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Creating the 3D head model from MRI... Please wait');

create_3d_model         % Creates a point cloud of the MRI scan

disp('Finished creating the final head model.')
disp(' ');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Point Cloud Registration - Iterative Closest Point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_model_and_scan     % Loads MRI model and scan point cloud
prepare_scan_for_icp    % Does preregistration
register_to_sphere      % Pre-registers onto a sphere
%%
user_rotate             % Makes user pre-align the scan
%filter_outside_sphere  % Filters points outside of a sphere around the scan point cloud
user_filter
register_to_mri         % Does final registraion onto MRI ptCloud
disp('Finished ICP registration.')
plot_aligned            % Plots registration results


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Registration of electrodes onto the registered scan model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_electrodes
load_cap
plot_electrodes
register_electrodes
project_electrodes
detect_electrodes
%extract_centroids





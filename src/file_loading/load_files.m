% GET PATH TO MRI SCAN - get full path to the mri scan
full_file_path = get_user_file_path('*.nii', 'Select the MRI NIfTI file');

% Load spm data from the path
header_info = spm_vol(full_file_path);
data_size = header_info.dim;
mri_t1 = spm_read_vols(header_info);
T_matrix = header_info.mat;  % transformation matrix
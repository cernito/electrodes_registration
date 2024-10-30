function [mri_t1, T_matrix] = load_files(full_file_path)
% Load spm data from the path

header_info = spm_vol(full_file_path);
data_size = header_info.dim;
mri_t1 = spm_read_vols(header_info);
T_matrix = header_info.mat;  % transformation matrix

end
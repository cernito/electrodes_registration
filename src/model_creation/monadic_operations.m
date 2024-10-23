
% Thresholding
th =30;
mri_t1_th = mri_t1 > th;
mid_slice_idx = size(mri_t1,3)/2;
mid_slice_th = mri_t1_th(:,:,mid_slice_idx);

% Get rid of components smaller than 100 pixels.
mri_t1_inverted = bwareaopen(imcomplement(mri_t1_th),10000,8);

mri_t1 = imcomplement(mri_t1_inverted);

% Fill holes
mri_t1_fill = imfill(mri_t1, 'holes');

% Dilatation
mri_t1 = bwareaopen(mri_t1_fill, 2000, 8);

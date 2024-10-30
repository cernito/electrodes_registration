function [pcHead_model, pcHead_scan] = load_model_and_scan()
    
    full_file_path = get_user_file_path('*.stl', 'Select the 3D head model made from MRI');
    mri_model = stlread(full_file_path);
    pcHead_model = pointCloud(mri_model.Points);
    
    
    full_file_path = get_user_file_path('*.stl', 'Select the scan of patients head');
    head_scan = stlread(full_file_path);
    pcHead_scan = pointCloud(head_scan.Points);

end
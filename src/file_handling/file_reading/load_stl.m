function head_model = load_stl()
    full_file_path = get_user_file_path('*.stl', 'Select the 3D head model made from MRI');

    if isequal(full_file_path, 0)
        head_model = [];
        return;
    end
    
    head_model = stlread(full_file_path);
end
function pcCap = load_cap()

file_path = get_user_file_path('.stl', 'cap model');
cap = stlread(file_path);
pcCap = pointCloud(cap.Points);

end
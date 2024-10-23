function full_file_path = get_user_file_path(file_type,message)
    
[file_name, file_path] = uigetfile(file_type, message);
if isequal(file_name,0) || isequal(file_path,0)
    disp('User canceled file selection.');
else
    disp(['User selected ', fullfile(file_path, file_name)]);
end
% Full path to the mri scan
full_file_path = fullfile(file_path, file_name);

end
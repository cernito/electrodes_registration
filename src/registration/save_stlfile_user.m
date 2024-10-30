function save_stlfile_user(file,file_name)
% Prompt the user for a file name and location
message = join(['Save ', file_name, ' as']);
[file_name, file_path] = uiputfile('*.stl', message);
if isequal(file_name,0) || isequal(file_path,0)
    disp('User canceled the save operation.');
else
    full_file_path = fullfile(file_path, file_name);
    
    % Now save the mesh as an STL file
    stlwrite(file, full_file_path);
    
    disp([join([file_name,' saved to ']), full_file_path]);
end
end
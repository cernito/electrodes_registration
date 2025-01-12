function full_file_path = get_user_file_path(file_type, message)
% Lets user pick a file of the given 'file_type' with the File Explorer.
%
% INPUT: 
%   file_type [string] ... '.stl', '.ply', '.elc', etc.
%   message [string] ... message displayed in the File Explorer
%
% OUTPUT:
%    full_file_path [string] ... string path to file

[file_name, file_path] = uigetfile(file_type, message);
if isequal(file_name,0) || isequal(file_path,0)
    disp('User canceled file selection.');
    full_file_path = 0;
    return
else
    disp(['User selected ', fullfile(file_path, file_name)]);
end

% Full path to the file 
full_file_path = fullfile(file_path, file_name);

end
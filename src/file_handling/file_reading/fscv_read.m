function coords = fscv_read(filename)
% READFCSV  Read a 3D Slicer .fcsv fiducial file and return x,y,z coords
%    coords = readFcsv(filename)
%
%    INPUT:
%        filename - path to the .fcsv file
%
%    OUTPUT:
%        coords   - N-by-3 matrix of [x, y, z] coordinates

    % Open the file in read-text mode
    fid = fopen(filename, 'rt');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    coords = [];  % Initialize an empty array for storing coordinates
    
    % Read file line-by-line
    while ~feof(fid)
        line = fgetl(fid);
        
        % Skip empty lines or lines starting with '#'
        if isempty(line) || startsWith(line, '#')
            continue
        end
        
        % Split the line by commas
        tokens = strsplit(line, ',');
        
        % Basic sanity check: ensure we have at least 5 columns
        % (the 1st is an ID, then x,y,z at indices 2,3,4)
        if numel(tokens) < 5
            continue
        end
        
        % Convert strings to numeric values
        x = str2double(tokens{2});
        y = str2double(tokens{3});
        z = str2double(tokens{4});
        
        % Append this row to coords
        coords = [coords; x, y, z]; %#ok<AGROW>
    end
    
    % Close the file
    fclose(fid);
end

%% Electrode detection

full_file_path = get_user_file_path('*.stl', 'Select the scan of patients head');
head_scan = stlread(full_file_path);
points = head_scan.Points;

epsilon = 5;
minPts = 10;
[idx, corePts] = dbscan(points, epsilon, minPts);
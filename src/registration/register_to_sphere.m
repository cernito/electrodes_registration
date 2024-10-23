% ICP parameters
metric = "pointToPoint";
extrapolate = true;
inlierRatio = 1;
maxIterations = 400;
tolerance = [0.001 0.01];

disp(''); disp('Performing the Point Cloud Registration of Head scan onto sphere... Please wait');

[tform1,alignedToSphere] = pcregistericp(pcHead_scan_filtered, pcSphere, ...
    Metric=metric, ...
    InlierRatio=inlierRatio, ...
    MaxIterations=maxIterations, ...
    Tolerance=tolerance, ...
    Extrapolate=extrapolate, ...
    Verbose=true);

disp('Finished first ICP registration.')

% fig = fig + 1;
% figure(fig); clf
% pcshowpair(alignedToSphere, pcMri_model)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% title('Registered scan to sphere')

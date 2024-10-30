function ptCloudAligned = register_to_mri(filteredCloud, pcMri_model)

    disp(''); disp('Performing the Point Cloud Registration of Head scan onto the MRI head model... Please wait');
    
    % ICP parameters
    metric = "pointToPoint";
    extrapolate = true;
    inlierRatio = 1;
    maxIterations = 400;
    tolerance = [0.001 0.01];
    
    
    % Downsampling
    %pcHead_scan = pcdownsample(pcHead_scan,"nonuniformGridSample",10);
    %pcMri_model = pcdownsample(pcMri_model,"nonuniformGridSample",10);
    
    % Performing ICP registrationpcshow(pcHead_scan.Location, labels)
    [tform2,ptCloudAligned] = pcregistericp(filteredCloud, pcMri_model, ...
        Metric=metric, ...
        InlierRatio=inlierRatio, ...
        MaxIterations=maxIterations, ...
        Tolerance=tolerance, ...
        Extrapolate=extrapolate, ...
        Verbose=true);

end
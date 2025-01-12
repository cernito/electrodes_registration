function pcFinal_positions = template_matching(scan, pcAligned_electrodes)

    % mainPath = fileparts(mfilename('fullpath'));
    % addpath(genpath(fullfile('C:\ČVUT\Bakalarka\electrodes_registration','src')));
    
    % ==========================
    % LOAD AND PREPROCESS THE HEAD SCAN
    % ==========================
    if isempty(scan)
        scan = load_stl
    end
    og_scan = scan;
    
    % ==========================
    % LOAD ALIGNED ELECTRODE POSITIONS
    % ==========================
    if isempty(pcAligned_electrodes)
        matFilePath = get_user_file_path('.mat',"Load file with electrode positions.");
        data = load(matFilePath,"pcElectrodes_aligned");
        pcAligned_electrodes = data.pcElectrodes_aligned; 
    end
    
    % =======================   =================================
    % SMOOTH AND REDUCE SCAN  |  ASSIGN CURVATURE VALUES TO POINT CLOUD COLOR
    % =======================   =================================
    
    scan = reduce_and_smooth(scan);
    %[GC_scan, MC_scan] = computeCurvaturesPara(scan, true);
    
    FV.faces = scan.ConnectivityList;
    FV.vertices = scan.Points;
    
    getderivatives = 0;
    [PrincipalCurvatures,~,~,~,~,~]= GetCurvatures( FV ,getderivatives);
    MC_scan = 0.5 * (PrincipalCurvatures(1,:) + PrincipalCurvatures(2,:))';
    
    plot_MC(scan, MC_scan);
    
    rgb_scan = MC_to_rgb(MC_scan);
    pcScan = pointCloud(scan.Points, "Color", rgb_scan);
    
    % ======================
    % ACCESS ALL TEMPLATES
    % ======================
    %%
    electrodeDir = uigetdir('..\data\templates\','Choose folder with template electrodes.');
    electrodeFiles = dir(fullfile(electrodeDir,'*.stl'));
    
    templatesAllPoints = [];
    facesList = cell(length(electrodeFiles),1);
    startIndices = zeros(length(electrodeFiles),1);
    endIndices = zeros(length(electrodeFiles),1);
    
    figure(1)
    hold on

    for i=1:length(electrodeFiles)
        electrodeFilePath = fullfile(electrodeDir, electrodeFiles(i).name);
        TRelectrode = stlread(electrodeFilePath);
    
        nPts = size(TRelectrode.Points,1);
        startIdx = size(templatesAllPoints,1) + 1;
        endIdx = startIdx + nPts - 1;

        templatesAllPoints = [templatesAllPoints; TRelectrode.Points];
        facesList{i} = TRelectrode.ConnectivityList;
    
        templates_centroids(i,:) = mean(TRelectrode.Points,1);
    
        startIndices(i) = startIdx;
        endIndices(i) = endIdx;
    end
    
    % =========================
    % ALLIGN TEMPLATES WITH FOUND ELECTRODE POSITIONS
    % =========================
    X = pcAligned_electrodes.Location;
    Y = templates_centroids;
    
    % Init full set of options %%%%%%%%%%
    opt.method='affine';    % use nonrigid registration
    
    opt.beta=2;             % the width of Gaussian kernel (smoothness)
    opt.lambda=2;           % regularization weight
    
    opt.viz=1;              % show every iteration
    opt.outliers = 0.1;     % noise weight
    opt.fgt=0;              % do not use FGT (default)
    opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
    opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)
    
    opt.max_it=100;         % max number of iterations
    opt.tol=1e-6;           % tolerance
    
    [Transform, C]=cpd_register(X,Y, opt);
    
    transformedPoints = cpd_transform(templatesAllPoints,Transform);
    
    % Rebuild each electrode as a triangulation
    templateElectrodes = cell(length(electrodeFiles),1);
    for i=1:length(electrodeFiles)
        pts = transformedPoints(startIndices(i):endIndices(i), :);
        faces = facesList{i};
        templateElectrodes{i} = triangulation(faces, pts);
    end
    
    %%
    % =====================
    % PREPROCEESS ORIGINAL SCAN
    % =====================
    redu_scan = reduce_tri(og_scan,0.99);
    pcRedu_scan = pointCloud(redu_scan.Points);
    
    %% TEMPLATE MATCHING LOOP
    
    figure(5); clf;
    trisurf(og_scan.ConnectivityList, og_scan.Points(:,1), og_scan.Points(:,2), og_scan.Points(:,3), 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
    axis equal;
    hold on;
    colormap gray;
    lighting gouraud; 
    camlight headlight;
    title('Head Scan with Matched Electrode Centroids');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    
    figure(6); clf;
    trisurf(og_scan.ConnectivityList, og_scan.Points(:,1), og_scan.Points(:,2), og_scan.Points(:,3), 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
    axis equal;
    hold on;
    colormap gray;
    lighting gouraud; 
    camlight headlight;
    title('Head Scan with Matched Electrode Centroids');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    
    % Ask user to choose computation option
    choice = questdlg('Do you want to use one template set or choose individual templates?', ...
        'Template Matching Options', ...
        'Use Whole Template Set', 'Choose Individual Templates', 'Use Whole Template Set');
    
    switch choice
        case 'Use Whole Template Set'
            template_matching_mode = 128;
        case 'Choose Individual Templates'
            template_matching_mode = 1;
    end

    choice = questdlg('Do you want to use parallel computation?', ...
        'Parallelization', ...
        'Yes', 'No', 'No');
    
    switch choice
        case 'Yes'
            parallel_compute = true;
        case 'No'
            parallel_compute = false;
    end

    if ~parallel_compute
        if template_matching_mode == 128
            % =============================
            % 128 TEMPLATES - LOOP OVER ALL ELECTRODE FILES
            % =============================
            
            results = zeros(length(electrodeFiles), 3);
            best_rmse = Inf(length(electrodeFiles), 1);
            
            availabilityMask = true(size(redu_scan.Points,1),1);
            exclusionRadius = 8;
        
            figure(9)
            hold on
    
            tic
            for i = 1:length(electrodeFiles)
                warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
                
                filename = split(electrodeFiles(i).name,'.');
                electrodeFilePath = fullfile(electrodeDir, filename(1));
                fprintf('Processing electrode %d\\%d: %s\n', i, length(electrodeFiles), filename{1});
                
                % Load and preprocess electrode
                electrode = templateElectrodes{i};
                electrode = reduce_and_smooth(electrode);
                
                % Map the MC values to RGB
                [GC_electrode, MC_electrode] = computeCurvature(electrode);
                % FV.faces = electrode.ConnectivityList;
                % FV.vertices = electrode.Points;
                % 
                % getderivatives = 0;
                % [PrincipalCurvatures,~,~,~,~,~]= GetCurvatures(FV,getderivatives);
                % MC_electrode = 0.5 * (PrincipalCurvatures(1,:) + PrincipalCurvatures(2,:))';
                % 
                rgb_electrode = MC_to_rgb(MC_electrode);
                pcElectrode = pointCloud(electrode.Points, "Color", rgb_electrode);
                
                % Move electrode to nearest point on reduced scan
                [pcMovedElectrode, centroid] = move_electrode_to_nearest_point(pcElectrode, redu_scan);
                
                % Segment part of the scan around the electrode
                radius = 20;
                scan_area = segment_scan_in_radius_with_mask(redu_scan, centroid, radius, availabilityMask);
                %scan_area = segment_scan_in_radius(redu_scan, centroid, radius);
                if isempty(scan_area) || size(scan_area.Points,1) < size(electrode.Points,1)
                    disp('No valid scan area found. Skipping this electrode.')
                    results(i,:) = pcAligned_electrodes.Location(i,:);
                    continue;
                end
                
                % Compute curvature values
                [GC_scan_area, MC_scan_area] = computeCurvature(scan_area);
                    
                % FV.faces = scan_area.ConnectivityList;
                % FV.vertices = scan_area.Points;
                % 
                % getderivatives = 0;
                % [PrincipalCurvatures,~,~,~,~,~]= GetCurvatures(FV,getderivatives);
                % MC_scan_area = 0.5 * (PrincipalCurvatures(1,:) + PrincipalCurvatures(2,:))';

                % Match scan_area curvature histogram to the template
                % 
                % MC_scan_area = match_hists(MC_scan_area, MC_electrode);
                
                % plot_MC(electrode, MC_electrode);
                % plot_MC(scan_area, MC_scan_area);

                % Map the MC values to RGB
                rgb_scan_area = MC_to_rgb(MC_scan_area);
                pcScan_area = pointCloud(scan_area.Points,"Color",rgb_scan_area);
                
                % Threshold values
                % [th_lower_mean, th_upper_mean] = get_mean_threshold(MC_electrode);
                % [th_lower_mode, th_upper_mode] = get_mode_threshold(MC_electrode);
                % th_lower = min(th_lower_mean, th_lower_mode);
                % th_upper = max(th_upper_mean, th_upper_mode);
                %[electrode, MC_electrode] = filter_scan_by_curvature(electrode, MC_electrode, 0.45, 0.55);
                
                % [th_lower_mean, th_upper_mean] = get_mean_threshold(MC_scan_area);
                % [th_lower_mode, th_upper_mode] = get_mode_threshold(MC_scan_area);
                % th_lower = min(th_lower_mean, th_lower_mode);
                % th_upper = max(th_upper_mean, th_upper_mode);
                %[scan_area, MC_scan_area] = filter_scan_by_curvature(scan_area, MC_scan_area, 0.45, 0.55);

                % Prealign electrode to the scan area plane
                pcAlignedElectrode = match_scan_plane(pcScan_area, pcMovedElectrode);
                % pcAlignedElectrode = pcMovedElectrode;

                min_val = min([MC_electrode; MC_scan_area],[],'omitnan');
                max_val = max([MC_electrode; MC_scan_area],[],'omitnan');
                rgb_electrode = MC_to_rgb_global(MC_electrode, min_val, max_val);
                rgb_scan_area = MC_to_rgb_global(MC_scan_area, min_val, max_val);
                pcAlignedElectrode.Color = rgb_electrode;
                pcScan_area.Color = rgb_scan_area;
                
                % Visualize the result
                % figure;
                % pcshow(pcScan_area);
                % title('Scan Area with Matched Curvature Values');
                % hold on
                % pcshow(pcAlignedElectrode);
                % hold off

                % Register electrode to scan area using ICP
                moving = pcAlignedElectrode;
                fixed = pcScan_area;
                
                moving.Normal = pcnormals(moving);
                fixed.Normal = pcnormals(fixed);
                
                [~,movingReg_orig,rmse_orig] = pcregistericp(moving,fixed,Metric="planeToPlaneWithColor", ...
                    MaxIterations=300, ...
                    Tolerance=[0.001 0.005]);
                
                moving = flip_electrode(moving);
                [~, movingReg_flipped, rmse_flipped] = pcregistericp(moving,fixed,Metric="planeToPlaneWithColor", ...
                    MaxIterations=300, ...
                    Tolerance=[0.001 0.005]);

                % Choose the best alignment
                if rmse_flipped < rmse_orig
                    rmse = rmse_flipped;
                    movingReg = movingReg_flipped;
                else
                    rmse = rmse_orig;
                    movingReg = movingReg_orig;
                end
                
                disp(['RMSE: ',num2str(rmse)]);
                disp(' ')
                
                % Show result
                % figure;
                % pcshow(movingReg, "MarkerSize", 600);
                % hold on;
                % pcshow(fixed);
                % title(['Aligned Electrode: ' electrodeFiles(i).name]);

                % final electrode centroid on the head
                centroid = mean(movingReg.Location, 1);
                
                % Plot the centroid on the same figure
                figure(5)
                hold on
                plot3(centroid(1), centroid(2), centroid(3), 'r.', 'MarkerSize', 30, 'LineWidth', 2);
                drawnow;

                % Save the aligned electrode points
                results(i,:) = centroid;
                best_rmse(i) = rmse;

                distances = vecnorm(redu_scan.Points - centroid, 2, 2);
                availabilityMask = availabilityMask & (distances > exclusionRadius);
    
            end
            toc
            warning('on','MATLAB:triangulation:PtsNotInTriWarnId');
        
            final_positions = results;
            figure(9)
            hold on
            pcshow(final_positions,'r','MarkerSize',600);
            pcshow(pcAligned_electrodes.Location,'b','MarkerSize',500)
            hold off
        
            figure(5)
            pcshow(pcAligned_electrodes.Location,'b','MarkerSize',500)
        
        end
        
        
        
        if template_matching_mode == 1
            % ==========================
            % 1 TEMPLATE FOR ALL ELECTRODES
            % ==========================
            
            % chosen_electrode = 23;
            
            % Load and preprocess electrode
            % electrodeFilePath = fullfile(electrodeDir, electrodeFiles(chosen_electrode).name);
            % electrode = templateElectrodes{chosen_electrode};
        
            electrode = load_stl
            electrode = reduce_and_smooth(electrode);
            %[GC_electrode, MC_electrode] = computeCurvature(electrode);
            FV.faces = electrode.ConnectivityList;
            FV.vertices = electrode.Points;
            
            getderivatives = 0;
            [PrincipalCurvatures,~,~,~,~,~]= GetCurvatures(FV,getderivatives);
            MC_electrode = 0.5 * (PrincipalCurvatures(1,:) + PrincipalCurvatures(2,:))';
            
            rgb_electrode = MC_to_rgb(MC_electrode);
            pcElectrode = pointCloud(electrode.Points, "Color", rgb_electrode);
            
            figure(99);  clf;
            
            availabilityMask = true(size(redu_scan.Points,1),1);
            exclusionRadius = 8;

            results = zeros(length(electrodeFiles), 3);
            best_rmse = Inf(length(electrodeFiles), 1);
            
            tic
            for i = 1:length(electrodeFiles)
                warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
                
                filename = split(electrodeFiles(i).name,'.');
                fprintf('Processing electrode %\\%d: %s\n',i,length(electrodeFiles), filename{1});
                %pcElectrode = move_template_to_aligned_electrode(pcElectrode, pcAligned_electrodes.Location(i,:));
                
                % Prealign electrode 
                electrodeForAlignment = templateElectrodes{i};
                pcElectrodeForAlignment = pointCloud(electrodeForAlignment.Points);
                pcElectrodeForAlignment.Normal = pcnormals(pcElectrodeForAlignment);
                pcElectrode.Normal = pcnormals(pcElectrode);
                [~, pcElectrode] = pcregistericp(pcElectrode,pcElectrodeForAlignment,"Metric","planeToPlane");
                
                % Move electrode to nearest point on reduced scan
                [pcMovedElectrode, centroid] = move_electrode_to_nearest_point(pcElectrode, redu_scan);
                %pcMovedElectrode = pcElectrode;
            
                % Segment part of the scan around the electrode
                radius = 22;
                scan_area = segment_scan_in_radius_with_mask(redu_scan, centroid, radius, availabilityMask);
                %scan_area = segment_scan_in_radius(redu_scan, centroid, radius);
                if isempty(scan_area) || size(scan_area.Points,1) < size(electrode.Points,1)
                    disp('No valid scan area found. Skipping this electrode.')
                    results(i,:) = pcAligned_electrodes.Location(i,:);
                    continue
                end
                
                % Compute the curvatrure
                %[GC_scan_area, MC_scan_area] = computeCurvature(scan_area);
                FV.faces = scan_area.ConnectivityList;
                FV.vertices = scan_area.Points;
                
                getderivatives = 0;
                [PrincipalCurvatures,~,~,~,~,~]= GetCurvatures(FV,getderivatives);
                MC_scan_area = 0.5 * (PrincipalCurvatures(1,:) + PrincipalCurvatures(2,:))';

                
                % [th_lower_mean, th_upper_mean] = get_mean_threshold(MC_scan_area);
                % [th_lower_mode, th_upper_mode] = get_mode_threshold(MC_scan_area);
                % th_lower = min(th_lower_mean, th_lower_mode);
                % th_upper = max(th_upper_mean, th_upper_mode);
                % [scan_area, MC_scan_area] = filter_scan_by_curvature(scan_area, MC_scan_area, th_lower, th_upper);
        

                % DO FOR MC VALUES
                % plot_MC(scan_area, MC_scan_area);
                rgb_scan_area = MC_to_rgb(MC_scan_area);
                pcScan_area = pointCloud(scan_area.Points,"Color",rgb_scan_area);
                
                % Prealign electrode to the scan area plane
                pcAlignedElectrode = match_scan_plane(pcScan_area, pcMovedElectrode);
                %pcAlignedElectrode = pcMovedElectrode;
                
                min_val = min([MC_electrode; MC_scan_area],[],'omitnan');
                max_val = max([MC_electrode; MC_scan_area],[],'omitnan');
                rgb_electrode = MC_to_rgb_global(MC_electrode, min_val, max_val);
                rgb_scan_area = MC_to_rgb_global(MC_scan_area, min_val, max_val);
                pcAlignedElectrode.Color = rgb_electrode;
                pcScan_area.Color = rgb_scan_area;

                % Register electrode to scan area using ICP
                moving = pcAlignedElectrode;
                fixed = pcScan_area;
                moving.Normal = pcnormals(moving,10);
                fixed.Normal = pcnormals(fixed,10);
                
                [~,movingReg_orig,rmse_orig] = pcregistericp(moving,fixed,Metric="planeToPlaneWithColor", ...
                    MaxIterations=300, ...
                    Tolerance=[0.001 0.005]);
                
                moving = flip_electrode(moving);
                [~, movingReg_flipped, rmse_flipped] = pcregistericp(moving,fixed,Metric="planeToPlaneWithColor", ...
                    MaxIterations=300, ...
                    Tolerance=[0.001 0.005]);
                
                % Choose the best alignment
                if rmse_flipped < rmse_orig
                    rmse = rmse_flipped;
                    movingReg = movingReg_flipped;
                else
                    rmse = rmse_orig;
                    movingReg = movingReg_orig;
                end
                
                disp(['RMSE: ',num2str(rmse)]);
                disp(' ')
                
                % Show result
                % figure(99);
                % pcshow(movingReg, "MarkerSize", 600);
                % hold on;
                % pcshow(fixed);
                % title(['Aligned Electrode: ' electrodeFiles(i).name]);
                % drawnow;
                
                % final electrode centroid on the head
                centroid = mean(movingReg.Location, 1);
                
                % Plot the centroid on the same figure
                figure(5)
                hold on
                plot3(centroid(1), centroid(2), centroid(3), 'r.', 'MarkerSize', 30, 'LineWidth', 2);
                drawnow;
            
                % Save the aligned electrode points
                results(i,:) = centroid;
                best_rmse(i) = rmse;

                distances = vecnorm(redu_scan.Points - centroid, 2, 2);
                availabilityMask = availabilityMask & (distances > exclusionRadius);
            
            end
            toc
            warning('on','MATLAB:triangulation:PtsNotInTriWarnId');
        end
        
        
        final_positions = results;
        figure(9)
        hold on
        pcshow(final_positions,'r','MarkerSize',600);
        pcshow(pcAligned_electrodes.Location,'b','MarkerSize',200)
        hold off
    
    %%
    else     
        % Disable warnings in curvature computation
        if template_matching_mode == 128
            figure(5); clf;
            trisurf(og_scan.ConnectivityList, og_scan.Points(:,1), og_scan.Points(:,2), og_scan.Points(:,3), 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
            axis equal;
            hold on;
            colormap gray;
            lighting gouraud; 
            camlight headlight;
            title('Head Scan with Matched Electrode Centroids');
            xlabel('X'); ylabel('Y'); zlabel('Z');

            results = zeros(length(electrodeFiles), 3);
            best_rmse = Inf(length(electrodeFiles), 1);
        
            parfor i = 1:length(electrodeFiles)
                warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
        
                electrodeFilePath = fullfile(electrodeDir, electrodeFiles(i).name);
                % fprintf('Processing electrode: %s\n', electrodeFiles(i).name);
            
                % Load and preprocess electrode
                electrode = templateElectrodes{i};
                electrode = reduce_and_smooth(electrode);
                [GC_electrode, MC_electrode] = computeCurvature(electrode);
                rgb_electrode = MC_to_rgb(MC_electrode);
                pcElectrode = pointCloud(electrode.Points, "Color", rgb_electrode);
            
                % Move electrode to nearest point on reduced scan
                [pcMovedElectrode, centroid] = move_electrode_to_nearest_point(pcElectrode, redu_scan);
            
                % Segment part of the scan around the electrode
                radius = 30;
                scan_area = segment_scan_in_radius(redu_scan, centroid, radius);
                if isempty(scan_area)
                    disp('No valid scan area found. Skipping this electrode.')
                end
            
                % Compute the curvatrure
                [GC_scan_area, MC_scan_area] = computeCurvature(scan_area);
                

                % DO FOR MC VALUES
                rgb_scan_area_MC = MC_to_rgb(MC_scan_area);
                pcScan_area_MC = pointCloud(scan_area.Points,"Color",rgb_scan_area_MC);

                % Prealign electrode to the scan area plane
                pcAlignedElectrode = match_scan_plane(pcScan_area_MC, pcMovedElectrode);
                
                
                min_val = min([MC_electrode; MC_scan_area],[],'omitnan');
                max_val = max([MC_electrode; MC_scan_area],[],'omitnan');
                rgb_electrode = MC_to_rgb_global(MC_electrode, min_val, max_val);
                rgb_scan_area = MC_to_rgb_global(MC_scan_area, min_val, max_val);
                pcAlignedElectrode.Color = rgb_electrode;
                pcScan_area_MC.Color = rgb_scan_area;

                % Register electrode to scan area using ICP
                moving = pcAlignedElectrode;
                fixed = pcScan_area_MC;
                moving.Normal = pcnormals(moving);
                fixed.Normal = pcnormals(fixed);
                
                [~,movingReg_orig,rmse_orig] = pcregistericp(moving,fixed,Metric="planeToPlaneWithColor", ...
                    MaxIterations=300, ...
                    Tolerance=[0.001 0.005]);
    
                moving = flip_electrode(moving);
                [~, movingReg_flipped, rmse_flipped] = pcregistericp(moving,fixed,Metric="planeToPlaneWithColor", ...
                    MaxIterations=300, ...
                    Tolerance=[0.001 0.005]);

                % Choose the best alignment
                if rmse_flipped < rmse_orig
                    rmse = rmse_flipped;
                    movingReg = movingReg_flipped;
                else
                    rmse = rmse_orig;
                    movingReg = movingReg_orig;
                end
        
                % final electrode centroid on the head
                centroid = mean(movingReg.Location, 1);
               
                results(i,:) = centroid;
                best_rmse(i) = rmse;

            end
            toc
            warning('on', 'curvatures:UnreferencedPoints');
                
            final_positions = results;
            figure(5)
            hold on
            pcshow(final_positions,'r','MarkerSize',600);
            hold off
            pcshow(pcAligned_electrodes.Location,'b','MarkerSize',500)

            
        elseif template_matching_mode == 1
            % ==========================
            % 1 TEMPLATE FOR ALL ELECTRODES
            % ==========================
            
            electrodes = {};
            while true
                % select an electrode STL file
                [filenames, pathname] = uigetfile({'*.stl'}, 'Select Electrode STL Files', 'MultiSelect', 'on');
        
                % Check if user canceled
                if isequal(filenames, 0)
                    break;
                end

                if ischar(filenames)
                    filenames = {filenames};
                end
                
                numElectrodes = length(filenames);
                for f_idx=1:numElectrodes
                    electrodeFilePath = fullfile(pathname, filenames{f_idx});
                    fprintf('Loading electrode: %s\n', electrodeFilePath);
                    
                    % Load and process electrode
                    electrode = stlread(electrodeFilePath); 
                    electrode = reduce_and_smooth(electrode);
                    [GC_electrode, MC_electrode] = computeCurvature(electrode);
                    rgb_electrode_MC = MC_to_rgb(MC_electrode);
                    rgb_electrode_GC = GC_to_rgb(GC_electrode);
                    pcElectrode = pointCloud(electrode.Points);
                    pcElectrode.Normal = pcnormals(pcElectrode);
                    electrodes{end+1} = {pcElectrode, rgb_electrode_MC, rgb_electrode_GC};
                end
            end
        
            if isempty(electrodes)
                disp('No electrodes were loaded.');
                pcFinal_positions = [];
                return;
            end
        
            numElectrodes = length(electrodes);
        
            
            results = zeros(length(electrodeFiles), 3);
            best_rmse = Inf(length(electrodeFiles), 1);
            
            electrodes_Const = parallel.pool.Constant(electrodes);

            tic
            for k = 1:numElectrodes
                figure(99);  clf;

                tic
                disp(['Processing with electrode ',num2str(k),'\',num2str(numElectrodes)])
                parfor i = 1:length(electrodeFiles)
                    warning('off','MATLAB:triangulation:PtsNotInTriWarnId');

                    % Prealign electrode 
                    electrodeForAlignment = templateElectrodes{i};
                    pcElectrodeForAlignment = pointCloud(electrodeForAlignment.Points);
                    pcElectrodeForAlignment.Normal = pcnormals(pcElectrodeForAlignment);
                    
                    % Load template
                    electrode_curr_arr = electrodes_Const.Value{k};
                    pcElectrode = electrode_curr_arr{1};
                    rgb_electrode_MC = electrode_curr_arr{2};
                    rgb_electrode_GC = electrode_curr_arr{3};
                    pcElectrode.Color = rgb_electrode_MC;
                    [~, pcElectrode] = pcregistericp(pcElectrode,pcElectrodeForAlignment,"Metric","planeToPlane");
                    
                    % Move electrode to nearest point on reduced scan
                    [pcMovedElectrode, centroid] = move_electrode_to_nearest_point(pcElectrode, redu_scan);
                    %pcMovedElectrode = pcElectrode;
                    
                    % Segment part of the scan around the electrode
                    radius = 30;
                    scan_area = segment_scan_in_radius(redu_scan, centroid, radius);
                    if isempty(scan_area) || size(scan_area.Points,1) < size(electrode.Points,1)
                        disp('No valid scan area found. Skipping this electrode.')
                        results(i,:) = pcAligned_electrodes.Location(i,:);
                        continue
                    end
                    
                    % Compute the curvatrure
                    [GC_scan_area, MC_scan_area] = computeCurvature(scan_area);
                    
                    % [th_lower_mean, th_upper_mean] = get_mean_threshold(MC_scan_area);
                    % [th_lower_mode, th_upper_mode] = get_mode_threshold(MC_scan_area);
                    % th_lower = min(th_lower_mean, th_lower_mode);
                    % th_upper = max(th_upper_mean, th_upper_mode);
                    % [scan_area, MC_scan_area] = filter_scan_by_curvature(scan_area, MC_scan_area, th_lower, th_upper);
                    

                    % DO FOR MC VALUES
                    % plot_MC(scan_area, MC_scan_area);
                    rgb_scan_area = MC_to_rgb(MC_scan_area);
                    pcScan_area = pointCloud(scan_area.Points,"Color",rgb_scan_area);
                    
                    % Prealign electrode to the scan area plane
                    pcAlignedElectrode = match_scan_plane(pcScan_area, pcMovedElectrode);
                    %pcAlignedElectrode = pcMovedElectrode;
                    
                    min_val = min([MC_electrode; MC_scan_area],[],'omitnan');
                    max_val = max([MC_electrode; MC_scan_area],[],'omitnan');
                    rgb_electrode = MC_to_rgb_global(MC_electrode, min_val, max_val);
                    rgb_scan_area = MC_to_rgb_global(MC_scan_area, min_val, max_val);
                    pcAlignedElectrode.Color = rgb_electrode;
                    pcScan_area.Color = rgb_scan_area;

                    % Register electrode to scan area using ICP
                    moving = pcAlignedElectrode;
                    fixed = pcScan_area;
                    moving.Normal = pcnormals(moving);
                    fixed.Normal = pcnormals(fixed);
                    
                    [~,movingReg_orig,rmse_orig] = pcregistericp(moving,fixed,Metric="planeToPlaneWithColor", ...
                        MaxIterations=300, ...
                        Tolerance=[0.001 0.005]);
        
                    moving = flip_electrode(moving);
                    [~, movingReg_flipped, rmse_flipped] = pcregistericp(moving,fixed,Metric="planeToPlaneWithColor", ...
                        MaxIterations=300, ...
                    Tolerance=[0.001 0.005]);
                    
                    % Choose the best alignment
                    if rmse_flipped < rmse_orig
                        rmse = rmse_flipped;
                        movingReg = movingReg_flipped;
                    else
                        rmse = rmse_orig;
                        movingReg = movingReg_orig;
                    end
                    
                    % Show result
                    % figure(99);
                    % pcshow(movingReg, "MarkerSize", 600);
                    % hold on;
                    % pcshow(fixed);
                    % title(['Aligned Electrode: ' electrodeFiles(i).name]);
                    % drawnow;
                    
                    % final electrode centroid on the head
                    centroid = mean(movingReg.Location, 1);
                
                    % Save the aligned electrode points
                    if rmse < best_rmse(i)
                        results(i,:) = centroid;
                        best_rmse(i) = rmse;
                    end

                end
                toc
                warning('on','MATLAB:triangulation:PtsNotInTriWarnId');
                figure(5); clf
                hold on
                plot3(results(:,1), results(:,2), results(:,3), 'r.', 'MarkerSize', 30, 'LineWidth', 2);
                drawnow;
            end

            final_positions = results;
        end
    end
    

    pcElectrodes_aligned = pointCloud(final_positions);
    save("aligned_electrodes.mat","pcElectrodes_aligned")
    

    % Review of RMSE scores
    figure(10); clf
    trisurf(og_scan.ConnectivityList, og_scan.Points(:,1), og_scan.Points(:,2), og_scan.Points(:,3), 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
    axis equal;
    hold on;
    good_rmse_mask = best_rmse < 0.5;
    plot3(results(good_rmse_mask,1), results(good_rmse_mask,2), results(good_rmse_mask,3), 'g.', 'MarkerSize', 30, 'LineWidth', 2);
    plot3(results(~good_rmse_mask,1), results(~good_rmse_mask,2), results(~good_rmse_mask,3), 'r.', 'MarkerSize', 30, 'LineWidth', 2);
    colormap gray;
    lighting gouraud; 
    camlight headlight;
    title('Results with rmse < 0.5');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    
    
    
    
    %% Final CPD align
    % elc_file = elc_read('template_elc.elc');
    % labels = elc_file.labels;
    % [~, sortOrder] = natsort(labels);
    % sorted_aligned_electrodes = pcAligned_electrodes.Location(sortOrder,:);

    X = final_positions(good_rmse_mask,:);
    Y = pcAligned_electrodes.Location; %sorted_aligned_electrodes(good_rmse_mask,:);

    figure;
    pcshow(X,'r','MarkerSize',600)
    hold on;
    pcshow(Y,'b','MarkerSize',600)
    
    % Init full set of options %%%%%%%%%%
    opt.method='nonrigid'; % use nonrigid registration
    
    opt.beta=6;            % the width of Gaussian kernel (smoothness)
    opt.lambda=15;          % regularization weight
    
    opt.viz=1;              % show every iteration
    opt.outliers = 1 - sum(good_rmse_mask) / length(good_rmse_mask)     % noise weight
    opt.fgt=0;              % do not use FGT (default)
    opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
    opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)
    
    opt.max_it=150;         % max number of iterations
    opt.tol=1e-10;           % tolerance
    
    figure(10)
    [Transform, ~]=cpd_register(X,Y, opt);
    
    pcAligned_electrodes = cpd_transform(pcAligned_electrodes.Location,Transform);
    
    figure(100); clf;
    pcshow(pcAligned_electrodes,'r','MarkerSize',600)
    hold on;
    %pcshow(scan.Points,[.7 .7 .7],'MarkerSize',200)
    trisurf(og_scan.ConnectivityList,og_scan.Points(:,1),og_scan.Points(:,2),og_scan.Points(:,3), ...
        'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');grid('on');
    axis('equal');
    camlight('headlight');
    lighting('gouraud');
    title('CPD alignment. Final positions.')
    
    
    disp(' ')
    disp('Finished template matching');


    pcFinal_positions = pointCloud(final_positions);
end

    










function T = nonlinear_contrast(vals)
    % Step 2: Apply non-linear contrast transformation
    alpha = 0.75;
    gamma = 1/(1-alpha);
    T = zeros(size(vals));

    half_idx = vals < 0.5;
    T(half_idx) = 0.5 * (2 * vals(half_idx)).^gamma;
    T(~half_idx) = 1 - 0.5 * (2 - 2 * vals(~half_idx)).^gamma;
end

function Mby3_rgb = MC_to_rgb_global(MC, min_val, max_val)
    % Step 1: Normalize curvature values to range [0, 1]
    MC(MC > 0.36) = 0.36;
    MC(MC < -0.36) = -0.36;
    min_val = -0.36;
    max_val = 0.36;
    normalized_MC = (MC - min_val) / (max_val - min_val);
    normalized_MC = max(0, min(1, normalized_MC)); % Clip values to [0, 1]
    
    %T = normalized_MC;
    T = nonlinear_contrast(normalized_MC);
    %high_val = T > 0.55;
    %low_val = T < 0.45;
    %T(high_val) = 1;
    %T(low_val) = 0;

    % Step 2: Map normalized values to the jet colormap
    num_colors = 256;                % Number of colors in colormap
    jet_colormap = turbo(num_colors);  % Generate jet colormap
    color_indices = round(T * (num_colors - 1)) + 1;
    rgb_values = jet_colormap(color_indices, :);
    
    % Step 3: Combine into M-by-3 RGB matrix
    Mby3_rgb = rgb_values;
end


function Mby3_rgb = MC_to_rgb(MC)
    % Given mean curvature values
    
    % Step 1: Normalize curvature values to range [0, 1]
    MC(MC > 0.36) = 0.36;
    MC(MC < -0.36) = -0.36;
    min_val = -0.36;  % Expected minimum
    max_val = 0.36;   % Expected maximum
    normalized_MC = (MC - min_val) / (max_val - min_val);
    normalized_MC = max(0, min(1, normalized_MC)); % Clip values to [0, 1]
    
    T = normalized_MC;
    %T = nonlinear_contrast(normalized_MC);
    % high_val = T > 0.55;
    % low_val = T < 0.45;
    % T(high_val) = 1;
    % T(low_val) = 0;

    % Step 2: Map normalized values to the jet colormap
    num_colors = 256;                % Number of colors in colormap
    turbo_colormap = turbo(num_colors);  % Generate jet colormap
    color_indices = round(T * (num_colors - 1)) + 1;
    rgb_values = turbo_colormap(color_indices, :);
    
    % Step 3: Combine into M-by-3 RGB matrix
    Mby3_rgb = rgb_values;
    
    % figure;
    % % Optional: Visualize the color values for the curvature
    % scatter3(1:length(MC), zeros(size(MC)), MC, ...
    %          50, Mby3_rgb, 'filled');
    % xlabel('Point Index'); ylabel('Y'); zlabel('Curvature');
    % title('Mean Curvature Colored by Jet Colormap');
    % colorbar;
end

function Mby3_rgb = GC_to_rgb(GC)
    % Step 1: Calculate adaptive min and max values (exclude outliers)
    min_val = prctile(GC, 1); % 1st percentile
    max_val = prctile(GC, 99); % 99th percentile
    
    % Step 2: Clamp GC values to the calculated range
    GC_clamped = max(min(GC, max_val), min_val);
    
    % Step 3: Normalize curvature values to range [0, 1]
    normalized_GC = (GC_clamped - min_val) / (max_val - min_val);
    normalized_GC = max(0, min(1, normalized_GC)); % Clip values to [0, 1]
    
    % Step 4: Map normalized values to the jet colormap
    num_colors = 256;                % Number of colors in colormap
    jet_colormap = jet(num_colors);  % Generate jet colormap
    color_indices = round(normalized_GC * (num_colors - 1)) + 1;
    rgb_values = jet_colormap(color_indices, :);
    
    % Step 5: Combine into M-by-3 RGB matrix
    Mby3_rgb = rgb_values;
    
    % figure;
    % % Optional: Visualize the color values for the curvature
    % scatter3(1:length(MC), zeros(size(MC)), MC, ...
    %          50, Mby3_rgb, 'filled');
    % xlabel('Point Index'); ylabel('Y'); zlabel('Curvature');
    % title('Mean Curvature Colored by Jet Colormap');
    % colorbar;
end

function [GC, MC] = computeCurvature(TR)
    [GC, MC] = curvatures(TR.Points(:, 1), TR.Points(:, 2), TR.Points(:, 3), TR.ConnectivityList);
end

function [GC, MC] = computeCurvaturesPara(scan_area, do_plot)
    % Setup parallel computation
    numChunks = 4; % Adjust based on available cores
    chunkIndices = round(linspace(1, size(scan_area.Points, 1) + 1, numChunks + 1)); % Chunk boundaries
    
    % Preallocate cell arrays for curvature results
    GC_chunks = cell(numChunks, 1);
    MC_chunks = cell(numChunks, 1);

    % Parallel computation
    parfor ic = 1:numChunks
        % Indices for the current chunk
        startIdx = chunkIndices(ic);
        endIdx = chunkIndices(ic + 1) - 1;
    
        % Extract vertices in the current chunk
        chunkVertices = scan_area.Points(startIdx:endIdx, :);
        chunkVertexIndices = startIdx:endIdx;
    
        % Filter connectivity list to include only triangles referencing chunk vertices
        faceMask = all(ismember(scan_area.ConnectivityList, chunkVertexIndices), 2); % Triangles fully within chunk
        filteredFaces = scan_area.ConnectivityList(faceMask, :);
    
        % Reindex the connectivity list to match chunk vertices
        [~, newIndices] = ismember(filteredFaces, chunkVertexIndices);
        filteredFaces = reshape(newIndices, size(filteredFaces));

        %validFaceMask = all(newIndices > 0, 2);
        %filteredFaces = newIndices(validFaceMask, :);

        if isempty(filteredFaces)
            GC_chunks{ic} = [];
            MC_chunks{ic} = [];
            continue
        end
        
        % Compute curvature for the current chunk
        [GC_chunk, MC_chunk] = curvatures(chunkVertices(:, 1), chunkVertices(:, 2), chunkVertices(:, 3), filteredFaces);
    
        % Store results in cell arrays
        GC_chunks{ic} = GC_chunk;
        MC_chunks{ic} = MC_chunk;
    end
    
    % Combine results from all chunks
    GC = vertcat(GC_chunks{:});
    MC = vertcat(MC_chunks{:});
    
    if do_plot
        % Plot
        tri = scan_area.ConnectivityList;
        x = scan_area.Points(:,1);
        y = scan_area.Points(:,2);
        z = scan_area.Points(:,3);
        
        img = figure(); clf;
        set(img, 'Position', [100 100 1200 600]); 
        
        subplot(1,2,1)
        hold on
        axis equal
        pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',GC,'FaceColor','interp','EdgeColor','none');
        clim([-0.1 0.1])
        colormap jet
        colorbar
        xlabel('x')
        ylabel('y')
        title('Estimated GC');
        hold off
    
        subplot(1,2,2)
        hold on
        axis equal
        pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',MC,'FaceColor','interp','EdgeColor','none');
        clim([-0.36 0.36])
        colormap jet
        colorbar
        xlabel('x')
        ylabel('y')
        title('Estimated MC');
        hold off
    end
end

function plot_MC(TR, MC)
    tri = TR.ConnectivityList;
    x = TR.Points(:,1);
    y = TR.Points(:,2);
    z = TR.Points(:,3);
    
    img = figure(); clf;
    set(img, 'Position', [100 100 1200 600]);
    
    hold on
    axis equal
    patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',MC,'FaceColor','interp','EdgeColor','none');
    clim([-0.5, 0.5]);
    colormap jet
    colorbar
    xlabel('x')
    ylabel('y')
    title('Estimated MC');
    hold off
end

function plot_GC(TR, MC)
    tri = TR.ConnectivityList;
    x = TR.Points(:,1);
    y = TR.Points(:,2);
    z = TR.Points(:,3);
    
    img = figure(); clf;
    set(img, 'Position', [100 100 1200 600]);
    
    hold on
    axis equal
    patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',MC,'FaceColor','interp','EdgeColor','none');
    clim([-0.1, 0.1]);
    colormap jet
    colorbar
    xlabel('x')
    ylabel('y')
    title('Estimated GC');
    hold off
end


function TR = reduce_tri(TR, reduceFactor)
    pTR.vertices = TR.Points;
    pTR.faces = TR.ConnectivityList;
    reduced_pTR = reducepatch(pTR,reduceFactor);
    new_TR.Points = reduced_pTR.vertices;
    new_TR.ConnectivityList = reduced_pTR.faces;
    TR = new_TR;
end


function TR = smooth_tri(TR, niter, strength)
    % Prepare the mesh structure for smoothing
    fvTR.vertices = TR.Points;
    fvTR.faces = TR.ConnectivityList;
    
    % Perform smoothing using the smoothpatch function
    smoothFVtri = smoothpatch(fvTR, 1, niter, strength);
    TR = triangulation(smoothFVtri.faces, smoothFVtri.vertices);
end


function TR = reduce_and_smooth(TR)
    TR = reduce_tri(TR,0.99);
    TR = smooth_tri(TR,1,0.2);
end


function pcAlignedElectrode = match_scan_plane(pcScan_area, pcElectrode)
    % =========================
    % PREALIGN TEMPLATE TO SCAN_AREA PLANE (WITH FLIP CORRECTION)
    % ==========================
   
    % Compute Scan Area Normal
    scanPoints = pcScan_area.Location;
    electrodePoints = pcElectrode.Location;
    electrodeCentroid = mean(electrodePoints, 1);
    % PCA to find local scan normal
    [coeff_scan, ~, ~] = pca(scanPoints);
    scanNormal = coeff_scan(:,3);

    % Ensure the scan normal points outward. If you have a known reference point, e.g., the head centroid:
    scanCentroid = mean(pcScan_area.Location,1);
    toCentroid = scanCentroid - electrodeCentroid;
    if dot(scanNormal, toCentroid) > 0
        % Flip to point outward (away from scan centroid)
        scanNormal = -scanNormal;
    end

    % Compute main (bump) direction of the electrode ---
    [coeff_elec, ~, ~] = pca(electrodePoints);
    electrodeBumpDirection = coeff_elec(:,3);

    % If the electrode bump direction doesn't face outward similarly, flip it
    if dot(electrodeBumpDirection, scanNormal) < 0
        electrodeBumpDirection = -electrodeBumpDirection;
    end

    % Rotate the electrode so that bump direction aligns with scan normal ---
    % Compute rotation axis and angle
    v = cross(electrodeBumpDirection, scanNormal); % Rotation axis
    s = norm(v);
    c = dot(electrodeBumpDirection, scanNormal);

    if s < 1e-9
        % Already aligned or opposite
        R = eye(3);
        if c < 0
            % Rotate by 180° if directly opposite
            R = rotation_matrix_about_axis([1;0;0], pi); % arbitrary axis
        end
    else
        v = v/s; % normalize rotation axis
        % Rodrigues' rotation formula
        K = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        R = eye(3) + K*sin(acos(c)) + K*K*(1 - c);
    end

    % Apply rotation about electrode centroid ---
    translated_points = electrodePoints - electrodeCentroid;
    rotated_points = (R * translated_points')';
    aligned_points = rotated_points + electrodeCentroid;

    pcAlignedElectrode = pointCloud(aligned_points, 'Color', pcElectrode.Color);
end


function flipped_electrode = flip_electrode(electrode)
    points = electrode.Location;
    centroid = mean(points, 1);

    translatedPoints = points - centroid;

    [coeff, ~, ~] = pca(translatedPoints);
    normalVector = coeff(:,3);

    % Reflection matrix about a plane with normal n:
    % R = I - 2*(n*n')
    normalVector = normalVector / norm(normalVector);
    reflectionMatrix = eye(3) - 2*(normalVector * normalVector');
    
    flippedPoints = (reflectionMatrix * translatedPoints')';

    finalPoints = flippedPoints + centroid;
    flipped_electrode = pointCloud(finalPoints, 'Color',electrode.Color,'Normal',electrode.Normal);
end


function R = rotation_matrix_about_axis(axis, theta)
    % ROTATION_MATRIX_ABOUT_AXIS: creates a rotation matrix for rotating by angle theta about a given axis.
    %
    % axis: [x; y; z] (should be normalized)
    % theta: rotation angle in radians
    %
    axis = axis / norm(axis);
    x = axis(1); y = axis(2); z = axis(3);
    c = cos(theta);
    s = sin(theta);
    R = [c+(1-c)*x*x,    (1-c)*x*y - s*z, (1-c)*x*z + s*y;
         (1-c)*y*x + s*z, c+(1-c)*y*y,    (1-c)*y*z - s*x;
         (1-c)*z*x - s*y, (1-c)*z*y + s*x, c+(1-c)*z*z];
end


function scan_area = segment_scan_in_radius(scan, centroid, radius)
    pcScan = pointCloud(scan.Points);
    [indices, ~] = findNeighborsInRadius(pcScan,centroid,radius);
    
    if isempty(indices) || size(indices,1) < 50
        disp('No neighbors found for electrode. Assigning template position.');
        scan_area = [];
        return
    end
    
    % Create a mask for filtering
    keepMask = false(size(scan.Points,1),1);
    keepMask(indices) = 1;
    
    % Filter the triangulation
    scan_area = filterTriangulation(scan, keepMask);
end

function scan_area = segment_scan_in_radius_with_mask(scan, centroid, radius, availabilityMask)
    % SEGMENT_SCAN_IN_RADIUS_WITH_MASK
    % Segments the scan within a given radius around a centroid, excluding points
    % marked in the exclusion mask.

    pcScan = pointCloud(scan.Points);

    % Find neighbors within the radius
    [indices, ~] = findNeighborsInRadius(pcScan, centroid, radius);

    % Exclude points that are within the exclusion mask
    validIndices = indices(availabilityMask(indices));

    if isempty(validIndices)
        scan_area = [];
        return;
    end

    % Create a mask for filtering
    keepMask = false(size(scan.Points, 1), 1);
    keepMask(validIndices) = true;

    % Filter the triangulation
    scan_area = filterTriangulation(scan, keepMask);
end




function [pcMovedElectrode, centroid] = move_electrode_to_nearest_point(pcElectrode, scan)
    templateElectrode = pcElectrode.Location;
    centroid = mean(templateElectrode,1);
    pcScan = pointCloud(scan.Points);
    
    % Step 1: Find the closest point on the scan to the centroid
    [closest_idx, ~] = findNearestNeighbors(pcScan, centroid, 1);
    closest_point = pcScan.Location(closest_idx, :);
    
    % Step 2: Compute the translation vector
    translation_vector = closest_point - centroid;
    
    % Step 3: Apply the translation to the electrode points
    movedElectrode = templateElectrode + translation_vector;
    
    % Step 4: Verify the new centroid matches the closest point
    centroid = mean(movedElectrode, 1);
    assert(norm(centroid - closest_point) < 1e-6, ...
        'Centroid alignment failed: Check the translation logic.');


    % Step 5: Update the electrode point cloud with the translated points
    pcMovedElectrode = pointCloud(movedElectrode, 'Color', pcElectrode.Color);
end

function [scan_filtered, MC_filtered] = filter_scan_by_curvature(scan, MC, th_lower, th_upper)
    % FILTER_SCAN_BY_CURVATURE
    % Removes points from the scan whose curvature is within [th_lower, th_upper].
    %
    % Inputs:
    %   scan      - triangulation object of the head scan
    %   MC        - mean curvature values (Nx1 array)
    %   th_lower  - lower threshold for curvature
    %   th_upper  - upper threshold for curvature
    %
    % Output:
    %   scan_filtered - filtered triangulation object
    
    % Create a logical mask: keep points outside the [th_lower, th_upper] range
    keepMask = (MC < th_lower) | (MC > th_upper);
    
    if sum(keepMask == 1) < 50
        [scan_filtered, MC_filtered] = filter_scan_by_curvature(scan, MC, th_lower + 0.02, th_upper - 0.02);
        return
    else
        % Filter the triangulation
        scan_filtered = filterTriangulation(scan, keepMask);
        MC_filtered = MC(keepMask);
    end
end

function [lower_mode, upper_mode] = get_mode_threshold(MC_scan)
    % Define number of bins for the histogram
    numBins = 100;
    
    % Compute histogram of curvature values
    [counts, edges] = histcounts(MC_scan, numBins);
    
    % Find the bin with the maximum count (mode)
    [~, maxBinIdx] = max(counts);
    
    % Compute the center of the mode bin
    th_mode = (edges(maxBinIdx) + edges(maxBinIdx + 1)) / 2;
    
    % Define a constant range around the mode (adjust 'const' as needed)
    const = 0.05;
    lower_mode = th_mode - const;
    upper_mode = th_mode + const;
    
    % Display the mode and threshold range
    % fprintf('Mode Curvature: %.4f\n', th_mode);
    % fprintf('Filtering Curvature Range: [%.4f, %.4f]\n', lower_mode, upper_mode);
    
    % Plot the histogram for visualization
    % figure;
    % histogram(MC_scan, numBins, 'FaceColor', 'b', 'EdgeColor', 'k');
    % hold on;
    % yLimits = ylim;
    % plot([lower_mode lower_mode], yLimits, 'r--', 'LineWidth', 2);
    % plot([upper_mode upper_mode], yLimits, 'r--', 'LineWidth', 2);
    % title('Histogram of Mean Curvature Values');
    % xlabel('Mean Curvature');
    % ylabel('Frequency');
    % legend('Curvature Distribution', 'Filtering Range');
    % hold off;
end

function [lower_mean, upper_mean] = get_mean_threshold(MC_scan)
    numBins = 100;

    % Compute mean curvature
    th_mean = mean(MC_scan,1,'omitnan');
    
    % Define a constant range around the mean (adjust 'const' as needed)
    const = 0.05; % Example value
    lower_mean = th_mean - const;
    upper_mean = th_mean + const;
    
    % Display the mean and threshold range
    % fprintf('Mean Curvature: %.4f\n', th_mean);
    % fprintf('Filtering Curvature Range: [%.4f, %.4f]\n', lower_mean, upper_mean);
    
    % Plot the histogram for visualization
    % figure;
    % histogram(MC_scan, numBins, 'FaceColor', 'b', 'EdgeColor', 'k');
    % hold on;
    % yLimits = ylim;
    % plot([lower_mean lower_mean], yLimits, 'g--', 'LineWidth', 2);
    % plot([upper_mean upper_mean], yLimits, 'g--', 'LineWidth', 2);
    % title('Histogram of Mean Curvature Values');
    % xlabel('Mean Curvature');
    % ylabel('Frequency');
    % legend('Curvature Distribution', 'Filtering Range');
    % hold off;
end




function pcMoved_template = move_template_to_aligned_electrode(pcTemplate, electrode_points)
    template_points = pcTemplate.Location;

    electrode_centroid = mean(electrode_points,1);
    template_centroid = mean(template_points, 1);

    translation_vector = electrode_centroid - template_centroid; 
    moved_template = template_points + translation_vector;
    pcMoved_template = pointCloud(moved_template,"Color",pcTemplate.Color);
end


function [sortedC, idx] = natsort(C)
    % Sort cell array lexicographically and numerically

    prefix = cell(size(C));
    numPart = zeros(size(C));
    
    % Loop through each string to extract parts
    for i = 1:length(C)
        tokens = regexp(C{i}, '([A-Za-z]+)(\d+)', 'tokens');
        if ~isempty(tokens)
            prefix{i} = tokens{1}{1};
            numPart(i) = str2double(tokens{1}{2});
        else
            % Handle cases without numbers if necessary
            prefix{i} = C{i};
            numPart(i) = 0;
        end
    end
    prefix_upper = upper(prefix);
    T = table(prefix_upper, numPart, 'VariableNames', {'Prefix', 'Number'});
    [sortedC, idx] = sortrows(T, {'Prefix', 'Number'}, {'ascend', 'ascend'});
end



function filteredTR = filterTriangulation(TR, keepMask)
    % disp('Filtering head scan triangulation...')
        
    % Ensure keepMask is logical and within bounds
    keepMask = logical(keepMask);
    if numel(keepMask) > size(TR.Points,1)
        keepMask = keepMask(1:size(TR.Points,1));
    elseif numel(keepMask) < size(TR.Points,1)
        % Pad keepMask with false for any points beyond the mask
        tmpMask = false(size(TR.Points,1),1);
        tmpMask(1:numel(keepMask)) = keepMask;
        keepMask = tmpMask;
    end
    
    % Indices of points to keep
    keepIndices = find(keepMask);
    
    % Filter the vertex list
    filteredVertices = TR.Points(keepIndices, :);

    % Filter faces: only keep faces whose all vertices are in keepIndices
    validFacesMask = all(ismember(TR.ConnectivityList, keepIndices), 2);
    filteredFaces = TR.ConnectivityList(validFacesMask, :);
    
    % Remap old vertex indices to new indices
    oldToNew = zeros(size(TR.Points,1),1);
    oldToNew(keepIndices) = 1:numel(keepIndices);
    filteredFaces = oldToNew(filteredFaces);

    % Create the new triangulation with updated faces and vertices
    filteredTR = triangulation(filteredFaces, filteredVertices);
    % disp('Completed triangulation filtering.')
    % disp(' ')
end
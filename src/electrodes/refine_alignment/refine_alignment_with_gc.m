function pcFinalElectrodePositions = refine_alignment_with_gc(scan, pcElectrodes_aligned)

    % ======================================================================
    % FIND PRECISE ELECTRODE POSITIONS ON HEAD SCAN USING CURVATURE FEATURES
    % ======================================================================
    %
    % DESCRIPTION:
    % Locates electrode positions by analysing curvature features
    % of a head scan triangulation.
    %
    % Segments the found high curvature clusters (which probably correspond
    % to electrode clusters) and uses aligned template to choose the 
    % fitting electrode position.
    %
    % INPUT:
    % - scan [triangulation]... STL file of the head scan 
    % - pcElectrodes_aligned [ptCloud] ... aligned template electrodes
    %
    % OUTPUT:
    % - pcFinalElectrodePositions [ptCloud] ... found electrode positions
    %
    %
    % PIPELINE:
    % 1. Performs curvature analysis on the head scan triangulation.
    %    For highly detailed scans it is better to smooth the mesh using
    %    taubin filter.
    % 2. Segments vertices with high curvature value. 
    %    By observation, high values correspond to electrode centers.
    % 3. [  ]
    % 4. For each template electrode find the closest point on the 
    %    segmented head scan point cloud.
    
    disp('Electrode position refinement.')
    disp('Finds precise electrode positions on head scan')
    disp('using curvature features.')
    
    fig = 1;
    close all;

    % =========================
    % 1. Add Required Paths
    % =========================

    mainPath = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(mainPath, '..', 'src'))); % Adjust as needed

    % Define paths for data and results folders
    dataPath = fullfile(mainPath,'..','..','..', 'data');
    resultsPath = fullfile(mainPath,'..','..','..', 'results');
    
    % Create results directory if it doesn't exist
    if ~exist(resultsPath, 'dir')
        mkdir(resultsPath);
    end

    % Explicitly add 'curvatures', 'smoothpatch_version1b', and 'CPD' folders
    %addpath(genpath(fullfile(mainPath, '..', 'curvatures')));
    %addpath(genpath(fullfile(mainPath, '..', 'smoothpatch_version1b')));
    %addpath(genpath(fullfile(mainPath, '..', 'CPD')));

    % !!!
    % ADD "curvatures" folder into MATLAB environment path.
    % ADD "smoothpatch_version1b" folder into MATLAB environment path.
    % ADD "CPD" folder into MATLAB environment path.
    % !!!
    
    % =========================
    % 2. Initialize Parallel Pool
    % =========================

    % Desired number of workers (adjust based on your system)

    % Check if a parallel pool is already running
    pool = gcp('nocreate'); % Gets the current pool without creating a new one

    if isempty(pool)
        try
            % Start a parallel pool with the specified number of workers
            pool = parpool('local');
            disp(['Parallel pool initialized with ', num2str(pool.NumWorkers), ' workers.']);
        catch ME
            warning(['Failed to initialize parallel pool: ', ME.message]);
            disp('Proceeding without parallel processing.');
            pool = [];
        end
    else
        disp(['Parallel pool already running with ', num2str(pool.NumWorkers), ' workers.']);
    end

    % =========================
    % 3. Load Scan and Electrode Positions
    % =========================

    if nargin < 1 || isempty(scan)
        scan = load_stl;
    end
    if nargin < 2 || isempty(pcElectrodes_aligned)
        matFilePath = get_user_file_path('.mat',"Load file with electrode positions.");
        data = load(matFilePath,"pcElectrodes_aligned");
        pcElectrodes_aligned = data.pcElectrodes_aligned;
        
        includes_nasion = false;
        if size(pcElectrodes_aligned.Location,1) == 129
            % !!! pcElectrodes_aligned also include nasion point on last
            % position !!!
            includes_nasion = true;
        end
    end

    % Load electrode labels
    templatePath = fullfile(dataPath, 'template_elc.elc');
    elc_electrodes = elc_read(templatePath);
    %elc_electrodes = elc_read(get_user_file_path('.elc','Load tempalte electrodes. Add the filepath to the template instead of this code.'));
    %elc_electrodes = elc_read("C:\ČVUT\Bakalarka\template_elc.elc");
    labels = elc_electrodes.labels;


    % =========================
    % 4. Plot the Current Scan for User Review
    % =========================

    figPlot = figure('Name', 'Review Scan', 'NumberTitle', 'off');
    trisurf(scan.ConnectivityList, scan.Points(:,1), scan.Points(:,2), scan.Points(:,3), ...
        'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');
    title('Current Head Scan');
    axis equal;
    grid on;
    camlight('headlight');
    lighting('gouraud');
    view(3);

    % =========================
    % 5. Segment the head cap from scan
    % =========================
    
    og_scan = scan;
    if includes_nasion
        pcElectrodes_aligned = pointCloud(pcElectrodes_aligned.Location(1:end-1,:));
        includes_nasion = false;
    end

    % dont include M1 and M2 and points around ears
    around_ear_idx = true(size(pcElectrodes_aligned.Location,1),1);
    around_ear_idx([13, 19, 65, 66, 105, 110]) = false;
    electrodes_without_m_electrodes = pcElectrodes_aligned.Location(around_ear_idx,:);
    pcElectrodes_without_m_electrodes = pointCloud(electrodes_without_m_electrodes);

    segmented_cap = segment_head_cap(scan, pcElectrodes_without_m_electrodes);
    scan = segmented_cap;

    % Remove nasion point
    %     if includes_nasion
    %         pcElectrodes_aligned = pointCloud(pcElectrodes_aligned.Location(1:end-1,:));
    %     end

    % =========================
    % 6. Reduce patch size to fasten curvature computation
    % =========================
    
    disp('Reducing patch size to fasten curvature computation... Please wait')
    
    pScan.vertices = scan.Points;
    pScan.faces = scan.ConnectivityList;
    
    % Number of faces (triangles) in the scan
    numFaces = size(scan.ConnectivityList, 1);
    
    
    % Set a dynamic reduction factor based on the number of vertices
    if numFaces < 250000
        reduceFactor = 0.9;      % Less reduction for smaller scans
    elseif numFaces < 450000
        reduceFactor = 0.6;     % Stronger reduction for detailed scans
    elseif numFaces < 600000
        reduceFactor = 0.4;     % Stronger reduction for high-resolution scans
    else 
        reduceFactor = 0.2;
    end
    
    tic
    reduced_pScan = reducepatch(pScan,reduceFactor);
    new_scan.Points = reduced_pScan.vertices;
    new_scan.ConnectivityList = reduced_pScan.faces;
    
    scan = new_scan;
    pcCap = pointCloud(scan.Points);
    
    disp('Finished patch reduction.'); disp(' ');
    toc

    tic;
    [GC_cap, MC_cap] = computeCurvatures(scan, true);
    toc
    fig = fig + 1;

    % =========================
    % 7. Ask User Whether to Smooth the Scan
    % =========================

    choice = questdlg('How do you want to smooth the scan before computing mean curvature?', ...
        'Smoothing Option', ...
        'Strong', 'Weak', 'None','None');

    did_smooth = true;
    switch choice
        case 'Strong'
            disp('User chose to smooth the scan.');
            disp('Smoothing the mesh to mitigate false high curvature values that are not electrodes...');
            
            % Prepare the mesh structure for smoothing
            FVscan.vertices = scan.Points;
            FVscan.faces = scan.ConnectivityList;
            
            % Perform smoothing using the smoothpatch function
            tic
            smoothFVscan = smoothpatch(FVscan, 1, 3);
            
            disp('Completed smoothing.');
            toc
            disp(' ');
            
            % Update the scan with the smoothed mesh
            scan = triangulation(smoothFVscan.faces, smoothFVscan.vertices);
            
            % Close the initial scan review figure and open a new one for the smoothed scan
            close(figPlot);
            figSmooth = figure('Name', 'Smoothed Scan', 'NumberTitle', 'off');
            trisurf(scan.ConnectivityList, scan.Points(:,1), scan.Points(:,2), scan.Points(:,3), ...
                'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');
            title('Smoothed Head Scan');
            axis equal;
            grid on;
            camlight('headlight');
            lighting('gouraud');
            view(3);

        case 'Weak'
            disp('User chose to smooth the scan.');
            disp('Smoothing the mesh to mitigate false high curvature values that are not electrodes...');
           
            % Prepare the mesh structure for smoothing
            FVscan.vertices = scan.Points;
            FVscan.faces = scan.ConnectivityList;
            
            tic
            smoothFVscan = smoothpatch(FVscan, 1, 1);

            disp('Completed smoothing.');
            toc
            disp(' ');

            % Update the scan with the smoothed mesh
            scan = triangulation(smoothFVscan.faces, smoothFVscan.vertices);

            % Close the initial scan review figure and open a new one for the smoothed scan
            close(figPlot);
            figSmooth = figure('Name', 'Smoothed Scan', 'NumberTitle', 'off');
            trisurf(scan.ConnectivityList, scan.Points(:,1), scan.Points(:,2), scan.Points(:,3), ...
                'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');
            title('Smoothed Head Scan');
            axis equal;
            grid on;
            camlight('headlight');
            lighting('gouraud');
            view(3);

        case 'None'
            disp('User chose not to smooth the scan.');
            disp('Proceeding with the original scan.');
            did_smooth = false;
            % Keep the original scan without smoothing
            % No action needed
    end
    
    % =========================
    % 8. Plot scan with electrode positions
    % =========================

    fig = fig + 1;
    figure(fig); clf
    trisurf(scan.ConnectivityList, scan.Points(:,1), scan.Points(:,2), scan.Points(:,3), ...
        'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');
    hold on
    pcshow(pcElectrodes_aligned.Location,'r','MarkerSize',600)
    title('Head Scan with initial electrode positions');
    axis equal;
    grid on;
    camlight('headlight');
    lighting('gouraud');
    view(3);
    
    % Ask user to choose computation option
    choice = questdlg('Do you want to locate electrodes on whole head at once or locate electrodes on prefound segments?', ...
        'Electrode Detection Option', ...
        'Whole head', 'By segments', 'Whole head');
    
    switch choice
        case 'Whole head'
            curvature_on_whole_head = true;
        case 'By segments'
            curvature_on_whole_head = false;
    end


    % =========================
    % 9. Plot Reduced and Filtered Scan
    % =========================

    fig = fig + 1; 
    figure(fig); clf;
    %pcshow(scan.Points)
    trisurf(scan.ConnectivityList,scan.Points(:,1),scan.Points(:,2),scan.Points(:,3),'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');
    grid('on');
    axis('equal');
    camlight('headlight');
    lighting('gouraud');

    % =========================
    % 10. Compute curvature features
    % =========================

    disp('Begining curvature computation... Please wait')
    
    globalTimerVal = tic;
    if curvature_on_whole_head
        % Compute curvatures
        
        if did_smooth
            tic
            [GC_cap, MC_cap] = computeCurvatures(scan, true);
            toc
        end

        % Compute curvatures paralellized
        fig = fig + 1;
        GC = GC_cap;
        MC = MC_cap;
        
        % =========================
        % 11. Segment points with high curvature value
        % =========================

        disp(' ');
        disp('Segmenting points with high curvature value.')
        disp('Finding the best segmentation parameters... Please wait')
    
        % Define a range of threshold values to test
        thresholdRange = 0.0:0.02:0.12; %linspace(prctile(MC, 75), prctile(MC, 95), 10);
        epsilonRange = 2:0.2:3.4; %2:1:8; % Range for DBSCAN epsilon <--- electrodes have a radius of circa 5mm
        expectedClusters = 128; % Target number of clusters
        maxClusterSize = 300; % Maximum allowable points in a valid cluster
        
        % Initialize variables to store results
        bestThreshold = 0;
        bestEpsilon = 0;
        bestNumClusters = 0;
        minClusterDifference = Inf;
        
        tic
        results = [];
        for threshold = thresholdRange
            for epsilon = epsilonRange
                % Segment points based on the current threshold
                selectedVertices = (GC < -threshold); % | (MC < -threshold);
                selectedPoints = scan.Points(selectedVertices, :);
                
                % Perform DBSCAN clustering
                minPoints = 3; % Minimum points per cluster
                [idx, ~] = dbscan(selectedPoints, epsilon, minPoints);
        
                % Count the valid clusters (exclude noise and oversized clusters)
                uniqueClusters = unique(idx);
                numValidClusters = 0;
                for clusterID = uniqueClusters'
                    if clusterID == -1 % Skip noise
                        continue;
                    end
                    
                    % Get points in cluster
                    clusterPoints = selectedPoints(idx == clusterID, :);
                    if size(clusterPoints,1) > maxClusterSize
                        continue
                    end

                    % Compute pairwise distances
                    pairwiseDistances = pdist(clusterPoints); % Pairwise Euclidean distances
                    
                    % Find the maximum distance (diameter of the cluster)
                    clusterDiameter = max(pairwiseDistances);
                    
                    % Discard clusters with large diameter size
                    if clusterDiameter < 10
                        otherClusterPoints = selectedPoints(idx ~= -1 & idx ~= clusterID, :);
                        if min(pdist2(clusterPoints, otherClusterPoints)) > 2
                            numValidClusters = numValidClusters + 1;
                        end
                    end
                    
                    % Count points in the current cluster
                    % clusterSize = sum(idx == clusterID); 
                    % if clusterSize <= maxClusterSize
                    %   numValidClusters = numValidClusters + 1;
                    % end
                end

                % Calculate the difference from the expected number of clusters
                clusterDifference = abs(numValidClusters - expectedClusters);
        
                % Store the results
                results = [results; threshold, epsilon, numValidClusters, clusterDifference];
        
                % Update the best parameters
                if clusterDifference <= minClusterDifference && numValidClusters < expectedClusters
                    minClusterDifference = clusterDifference;
                    bestThreshold = threshold;
                    bestEpsilon = epsilon;
                    bestNumClusters = numValidClusters;
                end
            end
        end
        toc
        disp('Finished parameter computation.')
    
        % Display the results
        disp(['Best Threshold: ', num2str(bestThreshold)]);
        disp(['Best Epsilon: ', num2str(bestEpsilon)]);
        disp(['Number of Valid Clusters: ', num2str(bestNumClusters)]);
        disp(' ')
        
        
        % Reshape results for visualization
        thresholds = unique(results(:, 1));
        epsilons = unique(results(:, 2));
        numClustersMatrix = reshape(results(:, 3), length(epsilons), length(thresholds));
        
        % Heatmap of number of clusters
        fig = fig + 1;
        figure(fig); clf;
        imagesc(thresholds, epsilons, numClustersMatrix);
        colorbar;
        xlabel('High MC Threshold');
        ylabel('Epsilon');
        title('Number of Clusters for Threshold and Epsilon');
        set(gca, 'YDir', 'normal'); % Correct orientation
        
        % ======================================================
        % Use the best threshold and epsilon
        % ======================================================
        disp('Using the best threshold value and epsilon for DBSCAN...')
    
        selectedVertices = (GC < -bestThreshold); %(MC < -bestThreshold);
        selectedPoints = scan.Points(selectedVertices, :);
        
        [idx, ~] = dbscan(selectedPoints, bestEpsilon, minPoints);
        
        % Filter valid clusters
        uniqueClusters = unique(idx);
        validClusters = {};
        validCentroids = {};
        for clusterID = uniqueClusters'
            if clusterID == -1 % Skip noise
                continue;
            end

            % Get points in the current cluster
            clusterPoints = selectedPoints(idx == clusterID, :);
            
            % Compute pairwise distances
            pairwiseDistances = pdist(clusterPoints); % Pairwise Euclidean distances
            
            % Find the maximum distance (diameter of the cluster)
            clusterDiameter = max(pairwiseDistances);
            
            % Discard clusters with large diameter size
            if clusterDiameter < 10
                otherClusterPoints = selectedPoints(idx ~= -1 & idx ~= clusterID, :);
                if min(pdist2(clusterPoints, otherClusterPoints)) > 2
                    validClusters{end+1} = clusterPoints;
                    validCentroids{end+1} = mean(clusterPoints,1);
                end
            end

            % if size(clusterPoints, 1) <= maxClusterSize
            %   validClusters{end+1} = clusterPoints;
            % end
        end
        
        numValidClusters = numel(validClusters);
        validPoints = vertcat(validCentroids{:});
        
        disp(['Number of Valid Clusters: ', num2str(numValidClusters)]);
        
        % Prepare colors
        clusterIDs = cell2mat(arrayfun(@(i) i * ones(size(validCentroids{i},1),1),1:numValidClusters, 'UniformOutput',false)');
        colors = lines(numValidClusters);
        pointColors = colors(clusterIDs,:);
        
        % Plot clusters
        fig = fig + 1; 
        figure(fig); clf;
        hold on;
        
        % Plot each valid cluster
        if numValidClusters > 0
           scatter3(validPoints(:, 1), validPoints(:, 2), validPoints(:, 3), ...
                    30, pointColors, 'filled');
        end
        
        % Plot selected points
        keepMask = false(size(scan.Points,1),1);
        keepMask(selectedVertices) = true;
        clusters = filterTriangulation(scan,keepMask);


        trisurf(clusters.ConnectivityList,clusters.Points(:,1),clusters.Points(:,2),clusters.Points(:,3),'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');
        tri = clusters.ConnectivityList;
        x = clusters.Points(:,1);
        y = clusters.Points(:,2);
        z = clusters.Points(:,3);
        pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',MC(selectedVertices),'FaceColor','interp','EdgeColor','none');
        caxis([-0.36 , 0.36])
        colormap jet
        colorbar

        % Add labels and finalize
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        title(['Filtered Clusters with Threshold = ', num2str(bestThreshold), ...
               ' and Epsilon = ', num2str(bestEpsilon)]);
        axis equal;
        grid on;
        hold off;
        disp(' ');
        
        disp('Completed.'); disp(' ')
        
        
        %% Find the closest cluster for each template electrode
        disp('Finding the closest cluster for each template electrode... Please wait')
        tic

        %elc_file = elc_read("C:\ČVUT\Bakalarka\template_elc.elc")
        %templateElectrodes = elc_file.pos;
        templateElectrodes = pcElectrodes_aligned.Location;
        if includes_nasion
            templateElectrodes = templateElectrodes(1:end-1,:);
        end
        segmentedPoints = validPoints;
        
        % Initialize array to store the nearest points
        nearestPoints = zeros(size(templateElectrodes)); % Nx3 for the nearest points
        nearestIndices = zeros(size(templateElectrodes, 1), 1); % Indices in segmentedPoints
        distances = zeros(size(templateElectrodes, 1), 1); % Distances to nearest points
        availabilityMask = true(size(segmentedPoints,1),1);

        % Create KD-Tree
        tree = KDTreeSearcher(segmentedPoints);

        averageSpacing = 20;
        thresholdDistance = averageSpacing * 0.8;
        exclusionRadius = floor(averageSpacing * 0.7);
        
        numElectrodes = size(templateElectrodes,1);
        
        % Sort Electrodes by closest distance and higher MC value
        minDistances = zeros(numElectrodes,1);
        maxCurvatureValue = zeros(numElectrodes,1);
        
        for i = 1:numElectrodes
            electrodePos = templateElectrodes(i, :);
            
            % Get closest distance
            % aboveZeroMask = segmentedPoints(:,3) > 0;
            dists = vecnorm(segmentedPoints - electrodePos, 2, 2);
            [minDistances(i), minIdx] = min(dists);
            
            % Get max MC value in small radius
            closestPoint = segmentedPoints(minIdx,:);
            smallRadius = 7;
            neighborsIdxCell = rangesearch(scan.Points, closestPoint, smallRadius);
            neighborsIdx = neighborsIdxCell{:};
            if isempty(neighborsIdx)
                maxCurvatureValue(i) = 0;
            else
                localMC = MC(neighborsIdx);
                localGC = GC(neighborsIdx);
                maxCurvatureValue(i) = max(localMC,[],'omitnan') + max(abs(localGC),[],'omitnan');
            end
        end
        
        % Min-Max normalization
        dMin = min(minDistances);
        dMax = max(minDistances);
        cMin = min(maxCurvatureValue);
        cMax = max(maxCurvatureValue);
        
        % Normalize
        dNorm = (minDistances - dMin) ./ (dMax - dMin);
        cNorm = (maxCurvatureValue - cMin) ./ (cMax - cMin);

        % Combine into score
        alpha = 1;
        beta = 1;
        score = alpha * (1 - dNorm) + beta * cNorm;

        % Sort electrodes by given score
        [~, sortedOrder] = sort(score, 'descend');
        
        sortedTemplateElectrodes = templateElectrodes(sortedOrder, :);
        sortedLabels = labels(sortedOrder);

        
        % ==========
        % Ensure the points are somewhat uniform across the head
        % ==========
        min_x = min(segmented_cap.Points(:,1));
        x_size = (max(segmented_cap.Points(:,1)) - min_x);
        leftMask = templateElectrodes(:,1) < (min_x + (0.5)*x_size);
        mid_x_mask = templateElectrodes(:,1) >= (min_x + (1/3)*x_size) & templateElectrodes(:,1) < (min_x + (2/3)*x_size);
        rightMask = templateElectrodes(:,1) > (min_x + (0.5)*x_size);
        
        min_y = min(segmented_cap.Points(:,2));
        y_size = (max(segmented_cap.Points(:,2)) - min_y);
        backMask = templateElectrodes(:,2) < (min_y + (1/3)*y_size);
        mid_y_mask = templateElectrodes(:,2) >= (min_y + (1/3)*y_size) & templateElectrodes(:,2) < (min_y + (2/3)*y_size);
        frontMask = templateElectrodes(:,2) > (min_y + (2/3)*y_size);

        top_z_threshold = min(segmented_cap.Points(:,3)) + (2/3) * max(segmented_cap.Points(:,3));
        topMask = templateElectrodes(:,3) > top_z_threshold;
        
        % For each region, find the sorted indices that belong to that region
        leftIdx  = find(leftMask(sortedOrder));
        rightIdx = find(rightMask(sortedOrder));
        topIdx   = find(topMask(sortedOrder));
        backIdx  = find(backMask(sortedOrder));
        midIdx = find(mid_x_mask(sortedOrder));
        
        % Minimum coverage
        minLeft  = 5;
        minRight = 5;
        minTop   = 5;
        minBack  = 5;
        minMid = 3;
        
        % Save found electrode positions (will be used for CPD registration)
        found_idx = false(size(templateElectrodes,1),1);
        
        % Select the top 18 electrodes
        selected_count = 0;
        
        % Search Scan Segments
        segmentTop.indices = topIdx
        segmentTop.minCount = minTop;
        segmentLeft.indices = leftIdx
        segmentLeft.minCount = minLeft;
        segmentRight.indices = rightIdx
        segmentRight.minCount = minRight;
        segmentBack.indices = backIdx
        segmentBack.minCount = minBack;
        segmentMid.indices = midIdx
        segmentMid.minCount = minMid;
        
        scan_segments = {
            segmentMid,...
            segmentLeft,...
            segmentRight,...
            segmentTop,...
            segmentBack
        };
        
        for idx = 1:length(scan_segments)
            positions_count = 0;

            segment = scan_segments{idx};
            found_indices = segment.indices;
            min_count = segment.minCount

            for i=1:length(found_indices)
                electrodePos = sortedTemplateElectrodes(found_indices(i), :);
                dists = vecnorm(segmentedPoints - electrodePos, 2, 2);
                withinThresholdIdx = find(dists <= thresholdDistance & availabilityMask);
                
                if isempty(withinThresholdIdx)
                    warning(['Electrode ', sortedLabels{found_indices(i)}, ' has no available points within threshold distance. Assigning to template position.']);
                    nearestPoints(found_indices(i), :) = [0,0,0]; %electrodePos;
                    distances(found_indices(i)) = Inf;
                    found_idx(found_indices(i)) = false;
                    continue;
                end
                
                %% Choose the point if it has high GC and MC value
                [minDist, relativeIdx] = min(dists(withinThresholdIdx));
                centroidIdx = withinThresholdIdx(relativeIdx);
                centroidPoint = segmentedPoints(centroidIdx, :);
    
                % Check MC values around centroid
                surroundedIdx = rangesearch(scan.Points, centroidPoint, 3);
                
                % Extract MC values of surrounded points
                candidate_MC = MC(surroundedIdx{:});
    
                % If no high value around centoid -> it is not electrode
                % position
                if max(candidate_MC) < 0.15 && min(candidate_MC) > -0.3
                    warning(['Electrode ', sortedLabels{found_indices(i)}, ' has no available points within threshold distance. Assigning to template position.']);
                    nearestPoints(found_indices(i), :) = [0,0,0]; %electrodePos;
                    distances(found_indices(i)) = Inf;
                    found_idx(found_indices(i)) = false;
                    continue
                end
            
                found_idx(found_indices(i)) = true;
                nearestPoints(found_indices(i), :) = centroidPoint;

                % Mark neighboring points within exclusion radius as unavailable
                [excludedIdx, ~] = rangesearch(tree, centroidPoint, exclusionRadius);
                availabilityMask(excludedIdx{:}) = false;
            
                fprintf('Electrode %s assigned to point %d with distance %.4f.\n', sortedLabels{i}, centroidIdx, minDist);
                
                selected_count = selected_count + 1;
                positions_count = positions_count + 1;
                if positions_count == min_count
                    break;
                end
                
            end
        end

        % Select remaining best positions
        for i = 1:size(sortedTemplateElectrodes,1)
            % Compute distances from the current template electrode to all
            % segmented cluster centroids
            electrodePos = sortedTemplateElectrodes(i, :);
            dists = vecnorm(segmentedPoints - electrodePos, 2, 2);
            
            % Find candidate points within threshold distance
            withinThresholdIdx = find(dists <= thresholdDistance & availabilityMask);
            
            if isempty(withinThresholdIdx)
                warning(['Electrode ', sortedLabels{i}, ' has no available points within threshold distance. Assigning to template position.']);
                nearestPoints(i, :) = [0,0,0]; %electrodePos;
                distances(i) = Inf;
                found_idx(i) = false;
                continue;
            end
            
            %% Choose the point if it has high GC and MC value
            [minDist, relativeIdx] = min(dists(withinThresholdIdx));
            centroidIdx = withinThresholdIdx(relativeIdx);
            centroidPoint = segmentedPoints(centroidIdx, :);

            % Check MC values around centroid
            surroundedIdx = rangesearch(scan.Points, centroidPoint,3);
            
            % Extract MC values of surrounded points
            candidate_MC = MC(surroundedIdx{:});

            % If no high value around centoid -> it is not electrode
            % position
            if max(candidate_MC) < 0.15 && min(candidate_MC) > -0.3
                warning(['Electrode ', sortedLabels{i}, ' has no available points within threshold distance. Assigning to template position.']);
                nearestPoints(i, :) = [0,0,0]; %electrodePos;
                distances(i) = Inf;
                found_idx(i) = false;
                continue
            end
            
            found_idx(i) = true;
            nearestPoints(i, :) = centroidPoint;
            
            % Mark neighboring points within exclusion radius as unavailable
            [excludedIdx, ~] = rangesearch(tree, centroidPoint, exclusionRadius);
            availabilityMask(excludedIdx{:}) = false;
            
            fprintf('Electrode %s assigned to point %d with distance %.4f.\n', sortedLabels{i}, centroidIdx, minDist);

            selected_count = selected_count + 1;
            if selected_count == 18
                break;
            end
        end

        

        
        % Restore original electrode order if necessary
        originalOrder = sortedOrder;
        unsortedNearestPoints = zeros(size(templateElectrodes));
        unsortedNearestPoints(originalOrder, :) = nearestPoints;
        unsortedDistances = Inf(size(templateElectrodes,1),1);
        unsortedDistances(originalOrder) = distances;
        found_idx(originalOrder) = found_idx;
        
        nearestPoints = unsortedNearestPoints;
        distances = unsortedDistances;


        disp('Finished nearest neighbor searching.'); 
        toc
        disp(' ')
        
        % Visualization
        fig = fig + 1; 
        figure(fig); clf;
        pcshow(nearestPoints,'r','MarkerSize',600)
        hold on;
        %pcshow(scan.Points,[.7 .7 .7],'MarkerSize',200)
        trisurf(scan.ConnectivityList,scan.Points(:,1),scan.Points(:,2),scan.Points(:,3),'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');
        grid('on');
        axis('equal');
        camlight('headlight');
        lighting('gouraud');
        
    end
    
    %% Compute curvatures on segmented part of head scan
    if (~curvature_on_whole_head)
        fig = fig + 1; 
        figure(fig);figure(fig+1);figure(fig+2);figure(fig+3);figure(fig+4);
        clf(fig); clf(fig+1); clf(fig+2); clf(fig+3); clf(fig+4)
        
        figure(fig+4)
        trisurf(og_scan.ConnectivityList,og_scan.Points(:,1),og_scan.Points(:,2),og_scan.Points(:,3), ...
            'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');
        hold on;
        grid('on');
        axis('equal');
        camlight('headlight');
        lighting('gouraud');


        % Get number of electrodes
        numElectrodes = size(pcElectrodes_aligned.Location,1);
        if includes_nasion
            numElectrodes = numElectrodes - 1;
        end
        

        % Segmentation area radius
        averageSpacing = 15;
        radius = ceil(averageSpacing * 5);
        exclusionRadius = ceil(averageSpacing * 0.4);
        
        % Initialize final electrode positions
        finalElectrodePositions = zeros(numElectrodes, 3);
        
        % Initialize cell array to store nearest points
        all_nearest_points = cell(numElectrodes, 1);

        % Initialize availability mask: tru indicated the point is
        % available
        availableMask = true(size(scan.Points, 1), 1);

        % Create a KD-Tree for efficient neighbor searching
        tree = KDTreeSearcher(scan.Points);

        for i=1:numElectrodes
            fprintf('Processing electrode %d/%d...\n', i, numElectrodes);

            % Template electrode position
            templateElectrode = pcElectrodes_aligned.Location(i, :);

            % Segment part of the scan around the electrode
            [indices, ~] = findNearestNeighbors(pcCap,templateElectrode,radius);

            if isempty(indices)
                disp(['No neighbors found for electrode ', num2str(i), '. Assigning template position.']);
                all_nearest_points{i} = templateElectrode;
                finalElectrodePositions(i, :) = templateElectrode;
                continue;
            end

            % Create a mask for filtering
            keepMask = false(size(scan.Points,1),1);
            keepMask(indices) = 1;

            % Filter the triangulation
            scan_area = filterTriangulation(scan, keepMask);


            % Check if scan_area is valid
            if isempty(scan_area.Points)
                disp(['No points found in scan area for electrode ', num2str(i), '. Assigning template position.']);
                all_nearest_points{i} = templateElectrode;
                finalElectrodePositions(i, :) = templateElectrode;
                continue;
            end
            
            % Compute curvatures
            [GC, MC] = computeCurvatures(scan_area, false);
            
            % ==========================
            % Selection of the point with maximal curvature
            % ==========================

            % Define curvature threshold
            highCurvatureThreshold = 0.08;

            % Identify points with curvature above the threshold
            highCurvatureMask = (MC > highCurvatureThreshold); %| (MC < -highCurvatureThreshold);
            highCurvaturePoints = scan_area.Points(highCurvatureMask, :);

            if isempty(highCurvaturePoints)
                disp(['No high curvature points found for Electrode ', num2str(i), '. Assigning template position.']);
                all_nearest_points{i} = templateElectrode;
                finalElectrodePositions(i, :) = templateElectrode;
                continue;
            end
            
            % Normalize curvature values (e.g., to range [0, 1])
            minCurvature = min(MC);
            maxCurvature = max(MC);
            normalizedCurvature = (MC - minCurvature) / (maxCurvature - minCurvature);
            
            % Compute distances to a template electrode
            distances = vecnorm(scan_area.Points - templateElectrode, 2, 2);
            
            % Normalize distances
            minDistance = min(distances);
            maxDistance = max(distances);
            normalizedDistance = (distances - minDistance) / (maxDistance - minDistance);
            
            % Extract normalized curvature and distance for selected vertices
            selectedCurvatures = normalizedCurvature(highCurvatureMask);
            selectedDistances = normalizedDistance(highCurvatureMask);            % Compute overall score

            % Compute score as weighted sum
            alpha = 0.8; % Weight for curvature
            beta = 0.2;  % Weight for distance
            scores = alpha * selectedCurvatures + beta * selectedDistances;

            % Sort points based on scores in descending order
            [~, sortedIdx] = sort(scores, 'descend');
            sortedHighCurvaturePoints = highCurvaturePoints(sortedIdx, :);
            sortedScores = scores(sortedIdx);
            
            % Initialize a flag to check if assignment is successful
            assignmentMade = false;
            
            % Iterate through sorted high-curvature points
            for j = 1:length(sortedScores)
                candidatePoint = sortedHighCurvaturePoints(j, :);
                
                % Find the index of the candidate point in the original scan
                [~, originalIdx] = ismember(candidatePoint, scan.Points, 'rows');
                
                if ~availableMask(originalIdx)
                    continue; % Skip if the point has been excluded
                end
                
                % Assign the candidate point to the current electrode
                all_nearest_points{i} = candidatePoint;
                finalElectrodePositions(i, :) = candidatePoint;
                assignmentMade = true;
                disp(' ')
                disp("FFOOOOOOOOOOFOOOOOOOOOOUUUUUUNNNDNNNNNDNDNNNNDNDDDDDDD")
                disp(' ')
                disp(' ')

                % Exclude all points within the exclusion radius of the assigned point
                [indicesWithinRadius, ~] = rangesearch(tree, candidatePoint, exclusionRadius);
                if ~isempty(indicesWithinRadius)
                    availableMask(indicesWithinRadius{:}) = false;
                end
                
                fprintf('Electrode %s assigned to point [%0.2f, %0.2f, %0.2f].\n', ...
                        labels{i}, candidatePoint(1), candidatePoint(2), candidatePoint(3));
                
                break; % Move to the next electrode after successful assignment
            end
            
            if ~assignmentMade
                % If no suitable high-curvature point was found, assign the template position
                disp(['No available high-curvature points for Electrode ', num2str(i), '. Assigning template position.']);
                all_nearest_points{i} = templateElectrode;
                finalElectrodePositions(i, :) = templateElectrode;
            end

            figure(fig+4);
            pcshow(finalElectrodePositions(i, :),'r','MarkerSize',600)
            hold on;

        end
    end

    
    %% TODO:
    % 1.1 Handle cases where more electrodes belong to same cluster.
    % 1.2 Or check if some electrode is not in expected position based on
    %     distances from neighboring electrodes and handle that.
    % 1.3 Or align again electrodes to these points and the electrodes which
    %     are not on positions update with the aligned position (SHOW which
    %     electrode were aligned like this)
    disp('Performing affine alignment of template electrodes to found nearestPoints.')
    disp('Reason: helps us find misaligned electrodes.')
    disp('Handling cases where some electrodes belong to same electrode on scan.')
    disp('If electrode is far from expected position, use its location from aligned template.')
    disp(' ')
    
    fig = 25;
    if ~curvature_on_whole_head
        nearestPoints = vertcat(all_nearest_points{:});
        
        found_idx = true(size(pcElectrodes_aligned.Location,1),1);
        if includes_nasion
            found_idx = found_idx(1:end-1);
        end
    end
    
    templateElectrodes = pcElectrodes_aligned.Location;
    if includes_nasion
        templateElectrodes = templateElectrodes(1:end-1,:);
    end

    % =======================
    % CPD Alignment on found electrodes
    % =======================
    
    X = nearestPoints(found_idx,:);
    
    whole_template = true;

    if whole_template
        around_ear_idx = true(size(pcElectrodes_aligned.Location,1),1);
        around_ear_idx([13, 19, 65, 66, 105, 110]) = false;
        templateElectrodesWithoutEars = templateElectrodes(around_ear_idx,:);
        Y = templateElectrodesWithoutEars;
    else
        around_ear_idx = found_idx;
        around_ear_idx([13, 19, 65, 66, 105, 110]) = false;
        templateElectrodesWithoutEars = templateElectrodes(around_ear_idx,:);
        Y = templateElectrodesWithoutEars;
    end
    
    % Init full set of options %%%%%%%%%%
    opt.method='nonrigid'; % use nonrigid registration
    
    opt.beta=2;            % the width of Gaussian kernel (smoothness)
    opt.lambda=10;          % regularization weight
    
    opt.viz=1;              % show every iteration
    if whole_template
        opt.outliers = 0.3; %1 - (sum(found_idx) / size(found_idx,1))    % noise weight
    else
        opt.outliers = 0.05
    end
    opt.fgt=0;              % do not use FGT (default)
    opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
    opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)
    
    opt.max_it=100;         % max number of iterations
    opt.tol=1e-10;           % tolerance
    
    [Transform, C]=cpd_register(X,Y, opt);

    alignedElectrodes = cpd_transform(pcElectrodes_aligned.Location,Transform);

    figure,cpd_plot_iter(X, Y); title('Before');
    figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');
    
    fig = fig + 1; 
    figure(fig); clf;
    pcshow(alignedElectrodes,'r','MarkerSize',600)
    hold on;
    %pcshow(scan.Points,[.7 .7 .7],'MarkerSize',200)
    trisurf(og_scan.ConnectivityList,og_scan.Points(:,1),og_scan.Points(:,2),og_scan.Points(:,3), ...
        'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');grid('on');
    axis('equal');
    camlight('headlight');
    lighting('gouraud');
    title('CPD')
    
    % =========================
    % Correct misaligned electrodes with aligned template electrode position
    % =========================

    disp('Correcting misaligned electrodes with aligned template electrode position...')
    tic

    % Parameters
    averageSpacing = 20; % Average Spacing of electrodes. 
                         % Also Maximum allowable distance from template 
                         % electrode to nearest point
    maxDistance = ceil(averageSpacing * 1);
    localRadius = averageSpacing * 0.2;
    
    % Initialize output array for final electrode positions
    finalElectrodePositions = nearestPoints; % Start with nearest points
    if includes_nasion
        finalElectrodePositions = [finalElectrodePositions; pcElectrodes_aligned.Location(end,:)];
    end
    
    % Add a label to misalligned electrodes for visualization
    missalignedIdx = false(size(templateElectrodes,1),1);
    betterMCpointIdx = false(size(templateElectrodes,1),1);
    
    % Check each template electrode
    for i = 1:size(templateElectrodes, 1)
        % Compute distance to the nearest point

        if all(nearestPoints(i, :) == [0,0,0])
            nearestPoints(i, :) = alignedElectrodes(i, :);
        end
        distanceToNearest = norm(nearestPoints(i, :) - alignedElectrodes(i, :));
        
        % If the nearest point is outside the diameter, replace it with the template position
        if distanceToNearest > maxDistance
            fprintf('Electrode %d: No nearest point within diameter. Using template position.\n', i);
            finalElectrodePositions(i, :) = alignedElectrodes(i, :);
            missalignedIdx(i) = true; 
            continue
        end

        % Find the closest point on scan
        dists = vecnorm(scan.Points - alignedElectrodes(i, :), 2, 2);
        closeIdx = dists <= maxDistance; % Points within allowable distance
        closePoints = scan.Points(closeIdx, :);
        closeMCValues = MC_cap(closeIdx); % Corresponding MC values
        closeGCValues = GC_cap(closeIdx);
        
        if isempty(closePoints)
            % Skip if no points are within the allowable distance
            fprintf('Electrode %d: No nearest point within diameter. Using template position.\n', i);
            finalElectrodePositions(i, :) = alignedElectrodes(i, :);
            missalignedIdx(i) = true;
            continue
        end
        
        % Find the MC value of the point in `nearestPoints(i, :)
        if found_idx(i)
            nearestPointIdx = find(ismember(scan.Points, nearestPoints(i, :), 'rows'), 1);

            if ~isempty(nearestPointIdx)
                % We found the exact match in the scan
                nearestPointMC = MC_cap(nearestPointIdx);
            else
                % The nearestPoints(i,:) wasn't found in scan.Points
                % -> Fallback to the minimal distance point
                [~, closestIdx] = min(dists);
                finalElectrodePositions(i, :) = scan.Points(closestIdx,:);
                nearestPointMC = MC_cap(closestIdx);
            end
        else
            [~, closestIdx] = min(dists);
            finalElectrodePositions(i, :) = scan.Points(closestIdx,:);
            nearestPointMC = MC_cap(closestIdx);
        end

        % Find the best MC value around the template position
        localDists = vecnorm(closePoints - alignedElectrodes(i, :), 2, 2); % Distances from template
        localIdx = localDists <= localRadius; % Indices of points within the local radius
        
        if any(localIdx)
            % Get the point with the highest MC value in the local radius
            [maxMCValue, maxLocalIdx] = max(closeMCValues(localIdx));
            [minMCValue, minLocalIdx] = min(closeMCValues(localIdx));
            [minGCValue, minGCLocalIdx] = min(closeGCValues(localIdx));
            bestMCPoints = closePoints(localIdx, :);
            bestGCPoints = closePoints(localIdx, :);
            maxMCPoint = bestMCPoints(maxLocalIdx, :); % Best point based on MC value
            minMCPoint = bestMCPoints(minLocalIdx, :);

            % Compare the best MC value around the template to the MC value of `nearestPoints(i, :)`
            % If the template point has higher MC value, choose it
            if (maxMCValue > nearestPointMC) && (minMCValue < -0.25)
                fprintf('Electrode %s: Template position has higher MC value. Using template position.\n', labels{i});
                if maxMCValue > nearestPointMC && maxMCValue > 0.3
                    finalElectrodePositions(i, :) = maxMCPoint; % Update to the point with the highest MC
                elseif minGCValue < -0.05
                    localGC = closeGValues(localIdx);
                    lowGCpointsIdx = localGC < -0.05;
                    minGCPoint = mean(bestGCPoints(lowGCpointsIdx,:), 1);
                    finalElectrodePositions(i, :) = minGCPoint; % Update to the point with the lowest MC
                elseif minMCValue < -0.25
                    %localMC = closeMCValues(localIdx);
                    %lowMCpointsIdx = localMC < -0.2;
                    %minMCPoint = mean(bestMCPoints(lowMCpointsIdx,:), 1);
                    finalElectrodePositions(i, :) = minMCPoint; % Update to the point with the lowest MC
                elseif maxMCValue > nearestPointMC
                    finalElectrodePositions(i, :) = maxMCPoint; % Update to the point with the highest MC
                end
                betterMCpointIdx(i) = true; % Mark as better point found
            else
                [~, closestIdx] = min(dists);
                finalElectrodePositions(i, :) = scan.Points(closestIdx,:);
            end
        end
    end

    disp('Completed correction.'); 
    toc
    disp(' ')
   

    % =====================
    % Save Aligned Electrodes
    % =====================
    pcElectrodes_aligned = pointCloud(finalElectrodePositions);
    % save("refined_electrodes.mat","pcElectrodes_aligned")
    % writeElectrodesToFcsv("..\results\final_electrode_positions.fcsv",labels,finalElectrodePositions);
    save(fullfile(resultsPath, 'refined_electrodes.mat'), 'pcElectrodes_aligned');
    writeElectrodesToFcsv(fullfile(resultsPath, 'final_electrode_positions.fcsv'), labels, finalElectrodePositions)


    visualize_correction()

    % =====================
    % Display the final electrode positions
    % =====================
    visualize_final_electrode_positions()

    pcFinalElectrodePositions = pointCloud(finalElectrodePositions);
    disp('Finished electrode refinement.'); 
    toc(globalTimerVal)
    disp(' ')



    % =================
    % Align Template Again
    % =================
    fig = fig + 1;
    figure(fig);

    X = finalElectrodePositions;
    Y = alignedElectrodes;
    if includes_nasion
        Y = Y(1:end-1,:);
    end
    
    % Init full set of options %%%%%%%%%%
    opt.method='nonrigid'; % use nonrigid registration
    
    opt.beta=4;            % the width of Gaussian kernel (smoothness)
    opt.lambda=10;          % regularization weight
    
    opt.viz=1;              % show every iteration
    opt.outliers = 0.12;       % noise weight
    opt.fgt=0;              % do not use FGT (default)
    opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
    opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)
    
    opt.max_it=100;         % max number of iterations
    opt.tol=1e-10;           % tolerance
    
    [Transform, C]=cpd_register(X,Y, opt);
    
    alignedElectrodes = cpd_transform(alignedElectrodes,Transform);

    figure,cpd_plot_iter(X, Y); title('Before');
    figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');
    
    fig = fig + 1; 
    figure(fig); clf;
    pcshow(alignedElectrodes,'r','MarkerSize',600)
    hold on;
    %pcshow(scan.Points,[.7 .7 .7],'MarkerSize',200)
    trisurf(og_scan.ConnectivityList,og_scan.Points(:,1),og_scan.Points(:,2),og_scan.Points(:,3), ...
        'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');grid('on');
    axis('equal');
    camlight('headlight');
    lighting('gouraud');
    title('CPD')
    

    % ====================
    % Find electrode (circle) shapes using RANSAC on prealigned positions
    % ====================
    
    % numElectrodes = size(alignedElectrodes,1);
    % otherElecMask = true(numElectrodes,1);
    % 
    % tree = KDTreeSearcher(scan.Points);
    % 
    % for i=1:numElectrodes
    % 
    %     % 1. Podívat se na okolí kolem elektrody - vybrat mesh.
    %     %   ? Poloměr okolí vybrat na základě vzdálenosti od ?
    %     %   ? nejbližší elektordy.                           ?
    %     %
    %     % 2. Vysegmentovat body s MC hodnotou nad daný threshold
    %     %   ? Použít pevný threshold nebo iterovat ?
    %     %
    %     % 3. V dané oblasti pomocí RANSAC najít kruh s podmínkou na poloměr 
    %     %    mezi 4-6 a cílem o kruh s největší průměrnou MC hodnotou.
    %     %
    %     % 4. Střed daného kruhu je nová pozice elektordy.
    %     %
    %     currentElectrode = alignedElectrodes(i,:);
    %     otherElecMask(i) = false;
    %     otherElectrodes = alignedElectrodes(otherElecMask,:);
    % 
    %     [closest_idx, distance] = knnsearch(otherElectrodes,currentElectrode,'K',1);
    %     searchRadius = distance * 0.9;
    % 
    %     [excludedIdx, ~] = rangesearch(tree, segmentedPoints(selectedIdx, :), exclusionRadius);
    %     availabilityMask(excludedIdx{:}) = false;
    % 
    %     otherElecMask(i) = true;
    % end
    % 
    % 




    
    % =======================================
    % NESTED FUNCTIONS
    % ======================================

    function [GC, MC] = computeCurvatures(scan_area, do_plot)
        % Setup parallel computation
        numChunks = pool.NumWorkers; % based on available cores
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
        GC = real(vertcat(GC_chunks{:}));
        MC = real(vertcat(MC_chunks{:}));
        
        if do_plot
            % Plot
            tri = scan_area.ConnectivityList;
            x = scan_area.Points(:,1);
            y = scan_area.Points(:,2);
            z = scan_area.Points(:,3);
            
            img = figure(fig); clf;
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
    
    function visualize_correction()
        % Visualization of updated positions
        fig = fig + 1; 
        figure(fig); clf;
        hold on;
        
        % Plot template electrodes
        h1 = scatter3(alignedElectrodes(:, 1), alignedElectrodes(:, 2), alignedElectrodes(:, 3), ...
                 100, 'g', 'filled', 'DisplayName', 'Template Electrodes');
        
        % Plot initial nearest points
        h2 = scatter3(nearestPoints(:, 1), nearestPoints(:, 2), nearestPoints(:, 3), ...
                 50, 'b', 'filled', 'DisplayName', 'Nearest Points');
        
        % Plot final positions
        h3 = scatter3(finalElectrodePositions(~missalignedIdx, 1), finalElectrodePositions(~missalignedIdx, 2), finalElectrodePositions(~missalignedIdx, 3), ...
                 100, 'r', 'filled', 'DisplayName', 'Final Electrode Positions');
        
        % Plot missaligned electrodes
        h4 = scatter3(finalElectrodePositions(missalignedIdx, 1), finalElectrodePositions(missalignedIdx, 2), finalElectrodePositions(missalignedIdx, 3), ...
                 100, 'y', 'filled', 'DisplayName', 'Missaligned Electrodes');
    
    
        % Add connections for replaced electrodes
        for i = 1:size(templateElectrodes, 1)
            if norm(nearestPoints(i, :) - templateElectrodes(i, :)) > maxDistance
                line([templateElectrodes(i, 1), finalElectrodePositions(i, 1)], ...
                     [templateElectrodes(i, 2), finalElectrodePositions(i, 2)], ...
                     [templateElectrodes(i, 3), finalElectrodePositions(i, 3)], 'Color', 'k', 'LineStyle', '--');
            end
        end
        
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        title('Electrode Positions: Template and Final');
        legend([h1, h2, h3, h4], 'Template Electrodes', 'Nearest Points', 'Final Electrode Positions', 'Missaligned Electrodes');
        axis equal;
        grid on;
        hold off;
    end

    function visualize_final_electrode_positions()
        fig = fig + 1; 
        figure(fig); clf;
        
        pcshow(finalElectrodePositions(~missalignedIdx & ~betterMCpointIdx,:),'g','MarkerSize',600)
        hold on;
        pcshow(finalElectrodePositions(missalignedIdx,:), 'r','MarkerSize',600)
        pcshow(finalElectrodePositions(betterMCpointIdx,:),[240 100 10]/256,'MarkerSize',600) 
        
        % Plot the head scan mesh
        trisurf(og_scan.ConnectivityList,og_scan.Points(:,1),og_scan.Points(:,2),og_scan.Points(:,3), ...
            'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');
        
        % Enhance lighting
        lighting('gouraud');
        camlight('left');
        camlight('right')
        
        % Other plot properties
        grid('on');
        axis('equal');
        title("Final detected electrodes")
        hold off;
    end

    function writeElectrodesToFcsv(filename, labels, points)
        % WRITEELECTRODESTOFCSV Write electrode positions to a Slicer Markups .fcsv file
        %
        %   filename - path to the .fcsv output file (e.g. "electrodes.fcsv")
        %   labels   - cell array of electrode labels (e.g. {'Fp1','Fpz','Fp2',...})
        %   points   - Nx3 matrix of final [X Y Z] electrode positions (in LPS)
        
        % Sanity checks
        if length(labels) ~= size(points,1)
            error('Number of labels (%d) must match number of electrode points (%d).',...
                length(labels), size(points,1));
        end
    
        % Attempt to open file for writing
        fid = fopen(filename, 'w');
        if fid == -1
            error('Could not open file "%s" for writing.', filename);
        end
    
        % Print the Slicer Markups header
        fprintf(fid, '# Markups fiducial file version = 5.4\n');
        fprintf(fid, '# CoordinateSystem = LPS\n');
        fprintf(fid, '# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n');
        
        % Each row: 
        %  id, x, y, z, ow, ox, oy, oz, vis, sel, lock, label, desc, associatedNodeID
        %
        % In example, each line ends with ",,,2,0" (thus, effectively 16 fields),
        % so we'll replicate that exact format (the extra columns aren't strictly
        % standard, but they mimic your sample output).
        %
        % e.g.
        % Fp1,42.2289,24.7369,316.4459,0,0,0,1,1,1,0,Fp1,,,2,0
        
        for i = 1:length(labels)
            thisLabel = labels{i};
            x = points(i,1);
            y = points(i,2);
            z = points(i,3);
            
            fprintf(fid, '%s,%.6f,%.6f,%.6f,0,0,0,1,1,1,0,%s,,,2,0\n', ...
                thisLabel, x, y, z, thisLabel);
        end
    
        fclose(fid);
        fprintf('Wrote %d electrodes to "%s"\n', length(labels), filename);
    end


end

function filteredTR = filterTriangulation(TR, keepMask)
    disp('Filtering head scan triangulation...')
        
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
    disp('Completed triangulation filtering.')
    disp(' ')
    
end

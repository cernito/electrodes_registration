function pcFinalElectrodePositions = refine_alignment_validation(scan, pcElectrodes_aligned)

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
    
    % close all;
    fig = 1;

    % =========================
    % 1. Add Required Paths
    % =========================

    mainPath = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(mainPath, '..', 'src'))); % Adjust as needed

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
    elc_electrodes = elc_read("C:\ČVUT\Bakalarka\template_elc.elc");
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
    segmented_cap = segment_head_cap(scan, pcElectrodes_aligned);
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
    if numFaces > 600000
        reduceFactor = 0.1;
    elseif numFaces > 500000
        reduceFactor = 0.3;     % Less reduction for smaller scans
    elseif numFaces > 400000
        reduceFactor = 0.5;     % Stronger reduction for high-resolution scans
    elseif numFaces > 300000
        reduceFactor = 0.6;
    else
        reduceFactor = 0.7;
    end
    
    reduced_pScan = reducepatch(pScan,reduceFactor);
    new_scan.Points = reduced_pScan.vertices;
    new_scan.ConnectivityList = reduced_pScan.faces;
    
    scan = new_scan;
    pcCap = pointCloud(scan.Points);
    
    disp('Finished patch reduction.'); disp(' ');
    
    tic
    [GC_cap, MC_cap] = computeCurvatures(scan, true);
    %     [GC, MC] = curvatures(scan.Points(:, 1), scan.Points(:, 2), scan.Points(:, 3), scan.ConnectivityList);
    %     toc
    % 
    %     tri = scan.ConnectivityList;
    %     x = scan.Points(:,1);
    %     y = scan.Points(:,2);
    %     z = scan.Points(:,3);
    %     
    %     img = figure(fig); clf;
    %     set(img, 'Position', [100 100 1200 600]); 
    %     
    %     subplot(1,2,1)
    %     hold on
    %     axis equal
    %     pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',GC,'FaceColor','interp','EdgeColor','none');
    %     caxis([-4 , -0.4])
    %     colormap jet
    %     colorbar
    %     xlabel('x')
    %     ylabel('y')
    %     title('Estimated GC');
    %     hold off
    % 
    %     subplot(1,2,2)
    %     hold on
    %     axis equal
    %     pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',MC,'FaceColor','interp','EdgeColor','none');
    %     caxis([-0.36 , 0.36])
    %     colormap jet
    %     colorbar
    %     xlabel('x')
    %     ylabel('y')
    %     title('Estimated MC');
    %     hold off



    fig = fig + 1;

    % =========================
    % 7. Ask User Whether to Smooth the Scan
    % =========================
    did_smooth = true;

    choice = questdlg('How do you want to smooth the scan before computing mean curvature?', ...
        'Smoothing Option', ...
        'Strong', 'Weak', 'None','None');

    switch choice
        case 'Strong'
            disp('User chose to smooth the scan.');
            disp('Smoothing the mesh to mitigate false high curvature values that are not electrodes...');
            
            % Prepare the mesh structure for smoothing
            FVscan.vertices = scan.Points;
            FVscan.faces = scan.ConnectivityList;
            
            % Perform smoothing using the smoothpatch function
            smoothFVscan = smoothpatch(FVscan, 1, 3);
            
            disp('Completed smoothing.');
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
            
            smoothFVscan = smoothpatch(FVscan, 1, 1);

            disp('Completed smoothing.');
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
    
    if curvature_on_whole_head
        % Compute curvatures
        tic
        [GC_cap, MC_cap] = computeCurvatures(scan, true);
        %[GC_cap, MC_cap] = curvatures(scan.Points(:, 1), scan.Points(:, 2), scan.Points(:, 3), scan.ConnectivityList);
        toc

        % Compute curvatures paralellized
        tic
        fig = fig + 1;
        GC = GC_cap;
        MC = MC_cap;
        toc
        
        % =========================
        % 11. Segment points with high curvature value
        % =========================

        disp(' ');
        disp('Segmenting points with high curvature value.')
        disp('Finding the best segmentation parameters... Please wait')
    
        

        % ================
        % LOW THRESHOLD SEARCH
        % ================

        % Define a range of threshold values to test
        thresholdRange = 0:0.01:0.1; % Adjust based on curvature values in your data
        epsilonRange = 2:0.2:4; %2:1:8; % Range for DBSCAN epsilon <--- electrodes have a radius of circa 5mm
        expectedClusters = 128; % Target number of clusters
        maxClusterSize = 30; % Maximum allowable points in a valid cluster
        
        % Initialize variables to store results
        bestThreshold = 0;
        bestEpsilon = 0;
        bestNumClusters = 0;
        minClusterDifference = Inf;
        
        results = [];
        for threshold = thresholdRange
            for epsilon = epsilonRange
                % Segment points based on the current threshold
                selectedVertices = (MC < -threshold); % | (MC < -threshold);
                selectedPoints = scan.Points(selectedVertices, :);
                
                % Perform DBSCAN clustering
                minPoints = 10; % Minimum points per cluster
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

                    % Compute pairwise distances
                    pairwiseDistances = pdist(clusterPoints,"fasteuclidean",CacheSize=100); % Pairwise Euclidean distances
                    
                    % Find the maximum distance (diameter of the cluster)
                    clusterDiameter = max(pairwiseDistances);

                    % Discard clusters with large diameter size
                    if clusterDiameter < 13
                        otherClusterPoints = selectedPoints(idx ~= -1 & idx ~= clusterID, :);
                        if min(pdist2(clusterPoints, otherClusterPoints,"fasteuclidean",CacheSize="maximal")) > 2
                            numValidClusters = numValidClusters + 1;
                        end
                    end

                    if numValidClusters > expectedClusters + 10
                        break;
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

        
        % ===================
        disp('Using the best threshold value and epsilon for DBSCAN...')
    
        selectedVertices = (MC < -bestThreshold); %(MC < -bestThreshold);
        selectedPoints = scan.Points(selectedVertices, :);
        
        [idx, ~] = dbscan(selectedPoints, bestEpsilon, minPoints);
        
        % Filter valid clusters
        uniqueClusters = unique(idx);
        validClustersLow = {};
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
            if clusterDiameter < 13
                otherClusterPoints = selectedPoints(idx ~= -1 & idx ~= clusterID, :);
                if min(pdist2(clusterPoints, otherClusterPoints)) > 2
                    validClustersLow{end+1} = clusterPoints;
                end
            end
            
            % if size(clusterPoints, 1) <= maxClusterSize
            %   validClustersLow{end+1} = clusterPoints;
            % end
        end
        
        numValidClusters = numel(validClustersLow);
        validPointsLow = vertcat(validClustersLow{:});
        disp(['Number of Valid Clusters: ', num2str(numValidClusters)]);
        
        % Prepare colors
        clusterIDs = cell2mat(arrayfun(@(i) i * ones(size(validClustersLow{i},1),1),1:numValidClusters, 'UniformOutput',false)');
        colors = lines(numValidClusters);
        pointColors = colors(clusterIDs,:);
        
        % Plot clusters
        fig = fig + 1;
        figure(fig); clf;
        hold on;
        
        % Plot each valid cluster
        if numValidClusters > 0
           scatter3(validPointsLow(:, 1), validPointsLow(:, 2), validPointsLow(:, 3), ...
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


        
        % ================
        % HIGH THRESHOLD SEARCH
        % ================

        % Define a range of threshold values to test
        thresholdRange = 0.02:0.02:0.3; % Adjust based on curvature values in your data
        epsilonRange = 2:0.2:4; %2:1:8; % Range for DBSCAN epsilon <--- electrodes have a radius of circa 5mm
        expectedClusters = 128; % Target number of clusters
        maxClusterSize = 30; % Maximum allowable points in a valid cluster
        
        % Initialize variables to store results
        bestThreshold = 0;
        bestEpsilon = 0;
        bestNumClusters = 0;
        minClusterDifference = Inf;

        results = [];
        for threshold = thresholdRange
            for epsilon = epsilonRange
                % Segment points based on the current threshold
                selectedVertices = (MC > threshold); % | (MC < -threshold);
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

                    % Compute pairwise distances
                    pairwiseDistances = pdist(clusterPoints); % Pairwise Euclidean distances
                    
                    % Find the maximum distance (diameter of the cluster)
                    clusterDiameter = max(pairwiseDistances);

                    % Discard clusters with large diameter size
                    if clusterDiameter < 3
                        otherClusterPoints = selectedPoints(idx ~= -1 & idx ~= clusterID, :);
                        if min(pdist2(clusterPoints, otherClusterPoints)) > 3
                            numValidClusters = numValidClusters + 1;
                        end
                    end
                    
                    if numValidClusters > expectedClusters + 10
                        break;
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
    
        selectedVertices = (MC > bestThreshold); %(MC < -bestThreshold);
        selectedPoints = scan.Points(selectedVertices, :);
        
        [idx, ~] = dbscan(selectedPoints, bestEpsilon, minPoints);
        
        % Filter valid clusters
        uniqueClusters = unique(idx);
        validClustersHigh = {};
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
            if clusterDiameter < 3
                otherClusterPoints = selectedPoints(idx ~= -1 & idx ~= clusterID, :);
                if min(pdist2(clusterPoints, otherClusterPoints)) > 3
                    validClustersHigh{end+1} = clusterPoints;
                end
            end

            % if size(clusterPoints, 1) <= maxClusterSize
            %   validClustersLow{end+1} = clusterPoints;
            % end
        end
        
        numValidClusters = numel(validClustersHigh);
        validPointsHigh = vertcat(validClustersHigh{:});
        disp(['Number of Valid Clusters: ', num2str(numValidClusters)]);
        
        % Prepare colors
        clusterIDs = cell2mat(arrayfun(@(i) i * ones(size(validClustersHigh{i},1),1),1:numValidClusters, 'UniformOutput',false)');
        colors = lines(numValidClusters);
        pointColors = colors(clusterIDs,:);
        
        % Plot clusters
        fig = fig + 1; 
        figure(fig); clf;
        hold on;
        
        % Plot each valid cluster
        if numValidClusters > 0
           scatter3(validPointsHigh(:, 1), validPointsHigh(:, 2), validPointsHigh(:, 3), ...
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

        
        % ==========================
        % VALIDATION OF CLUSTERS
        % ==========================
        
        disp('Validating high-curvature clusters with low-curvature clusters...');
        validClusters = {};
        
        % Validate clusters
        lowClusterCentroids = cellfun(@mean, validClustersLow, 'UniformOutput', false);
        lowClusterCentroids = vertcat(lowClusterCentroids{:});
        kdTreeLow = KDTreeSearcher(lowClusterCentroids);
        validationRadius = 4;
        
        for i = 1:numel(validClustersHigh)
            clusterHigh = validClustersHigh{i};
        
            % Compute the centroid of the high-curvature cluster
            centroidHigh = mean(clusterHigh, 1);
        
            % Query KD-tree for the nearest low-curvature cluster
            [~, distLowNearest] = knnsearch(kdTreeLow, centroidHigh);
        
            % Check if a nearby low-curvature cluster exists within the validation radius
            if distLowNearest <= validationRadius
                % Cluster is valid
                validClusters{end+1} = clusterHigh;
            end
        end
        
        numValidClusters = numel(validClusters);
        disp(['Number of Valid Clusters: ', num2str(numel(validClusters))]);
        
        % Combine valid clusters into a single set of points
        if ~isempty(validClusters) && numValidClusters > 17
            validPoints = vertcat(validClusters{:});
            
            % Map valid points back to the original vertex indices
            validIndices = ismember(scan.Points, validPoints, 'rows');
    
            % Update the selectedVertices mask
            selectedVertices = validIndices;
        else
            validPoints = validPointsHigh;
            validClusters = validClustersHigh;
            numValidClusters = numel(validClusters);
            disp('No validated combined-curvature clusters found.');
        end


        % Prepare colors
        clusterIDs = cell2mat(arrayfun(@(i) i * ones(size(validClusters{i},1),1),1:numValidClusters, 'UniformOutput',false)');
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

        
        %% Plot aligned electrodes
        pc_elec = pcElectrodes_aligned;
        [tform,aligned,rmse]=pcregistericp(pc_elec,pointCloud(selectedPoints))
        
        fig = fig + 1;
        figure(fig); clf;
        pcshow(validPoints,'r',"MarkerSize",400)
        hold on
        
        pcshow(pc_elec.Location,"g","MarkerSize",80)
        
        pcshow(aligned.Location,"b","MarkerSize",100)
        hold off
        
        %% Find the closest cluster for each template electrode
        disp('Finding the closest cluster for each template electrode... Please wait')
    
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
        thresholdDistance = averageSpacing * 2.3;
        exclusionRadius = floor(averageSpacing * 0.7);

        % Calculate Available Points per Electrode
        nearestNeighborCounts = zeros(size(templateElectrodes,1),1);
        for i = 1:size(templateElectrodes,1)
            electrodePos = templateElectrodes(i, :);
            dists = vecnorm(segmentedPoints - electrodePos, 2, 2);
            nearestNeighborCounts(i) = sum(dists <= averageSpacing);
        end
        %nearestNeighborCounts(nearestNeighborCounts == 0) = inf;

        % Sort electrodes by Increasing Number of Available Points
        [~, sortedOrder] = sort(nearestNeighborCounts, 'descend');
        sortedTemplateElectrodes = templateElectrodes(sortedOrder, :);
        sortedLabels = labels(sortedOrder);

        % Save found electrode positions (will be used for CPD registration)
        found_idx = true(size(templateElectrodes,1),1);

        for i = 1:size(templateElectrodes,1)
            % Compute distances from the current template electrode to all segmented points
            electrodePos = templateElectrodes(i, :);
            dists = vecnorm(segmentedPoints - electrodePos, 2, 2);
            
            % Find candidate points within threshold distance
            withinThresholdIdx = find(availabilityMask);
            
            if isempty(withinThresholdIdx)
                warning(['Electrode ', sortedLabels{i}, ' has no available points within threshold distance. Assigning to template position.']);
                nearestPoints(i, :) = [0,0,0]; %electrodePos;
                distances(i) = Inf;
                found_idx(i) = false;
                continue;
            end
            
            %% Find the nearest available point
            [minDist, relativeIdx] = min(dists(withinThresholdIdx));
            %selectedIdx = withinThresholdIdx(relativeIdx);
            
            % Find the point with highest MC value
            %[maxValue, relativeIdx] = max(MC(withinThresholdIdx));
            %selectedIdx = withinThresholdIdx(relativeIdx);
            
            % Extract MC values for candidate points
            candidate_MC = MC(withinThresholdIdx);

            % Extract distances for candidate points
            candidate_dist = dists(withinThresholdIdx);

            % Normalize distance and MC values (min-max normalization)
            normalized_dist = (candidate_dist - min(candidate_dist)) / (max(candidate_dist) - min(candidate_dist));
            normalized_dist = 1 - normalized_dist;
            normalized_MC = (candidate_MC - min(candidate_MC)) / (max(candidate_MC) - min(candidate_MC));

            % Handle cases where all distances or MC values are the same
            if all(normalized_dist == normalized_dist(1))
                normalized_dist(:) = 0.5;  % Neutral value
            end
            
            if all(normalized_MC == normalized_MC(1))
                normalized_MC(:) = 0.5;    % Neutral value
            end

            % Compute combined score
            weight_dist = 0.2;
            weight_MC = 0.8;
            combined_score = (weight_dist * normalized_dist) + (weight_MC * normalized_MC);

            % Select the point with the highest combined score
            [maxScore, relativeIdx] = max(combined_score);
            selectedIdx = withinThresholdIdx(relativeIdx);

            % Store the nearest point, index, and distance
            nearestPoints(i, :) = segmentedPoints(selectedIdx, :);
            nearestIndices(i) = selectedIdx;
            distances(i) = minDist;

            % Mark neighboring points within exclusion radius as unavailable
            [excludedIdx, ~] = rangesearch(tree, segmentedPoints(selectedIdx, :), exclusionRadius);
            availabilityMask(excludedIdx{:}) = false;

            fprintf('Electrode %s assigned to point %d with distance %.4f.\n', sortedLabels{i}, selectedIdx, minDist);
        end
        
        % Restore original electrode order if necessary
        % unsortedNearestPoints = zeros(size(templateElectrodes));
        % unsortedNearestPoints(sortedOrder, :) = nearestPoints;
        % unsortedDistances = Inf(size(templateElectrodes,1),1);
        % unsortedDistances(sortedOrder) = distances;
        % 
        % sortedFoundIdx = found_idx;
        % inverseOrder = zeros(size(sortedOrder));
        % inverseOrder(sortedOrder) = 1:length(sortedOrder);
        % found_idx = sortedFoundIdx(inverseOrder);
        % %found_idx(originalOrder) = found_idx;
        % 
        % nearestPoints = unsortedNearestPoints;
        % distances = unsortedDistances;

        disp('Finished nearest neighbor searching.'); disp(' ')
        
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


            %             % Define parameter ranges for clustering
            %             thresholdRange = 0.02:0.02:0.3; % Adjust based on curvature values in your data
            %             epsilonRange = 2:0.2:5; % Range for DBSCAN epsilon
            %             maxClusterSize = 10; % Maximum allowable points in a valid cluster
            %             expectedClusters = 1;
            %             minPoints = 3; % Minimum points per cluster for DBSCAN
            % 
            %             % Generate all combinations using combvec and transpose
            %             param_combinations = combvec(thresholdRange, epsilonRange)';
            %             num_combinations = size(param_combinations, 1);
            %             
            %             % Initialize score vector
            %             scores = -Inf(num_combinations,1);
            %             
            %             % Precompute max curvature and max distance for normalization
            %             maxCurvature = max(abs(MC));
            %             maxDistance = max(pdist2(scan_area.Points, scan_area.Points, 'euclidean'), [], 'all');
            % 
            %             % Define weights for curvature and distance
            %             weight_curvature = 0.8;
            %             weight_distance = 0.2;
            %             
            %             parfor k = 1:num_combinations
            %                 % Extract current threshold and epsilon
            %                 threshold = param_combinations(k, 1)
            %                 epsilon = param_combinations(k, 2)
            % 
            %                 % Segment points based on the current threshold
            %                 selectedVertices = (MC > threshold); %| (MC < -threshold);
            %                 selectedPoints = scan_area.Points(selectedVertices, :);
            %                 selectedCurvatures = MC(selectedVertices);
            %                 
            %                 if isempty(selectedPoints)
            %                     disp('No points segmented with threshold.')
            %                     continue
            %                 end
            %         
            %                 % Perform DBSCAN clustering
            %                 [idx, ~] = dbscan(selectedPoints, epsilon, minPoints);
            %         
            %                 % Count the valid clusters (exclude noise and oversized clusters)
            %                 uniqueClusters = unique(idx);
            %                 numValidClusters = 0;
            %                 for clusterID = uniqueClusters'
            %                     if clusterID == -1 % Skip noise
            %                         continue;
            %                     end
            %                     % Count points in the current cluster
            %                     clusterSize = sum(idx == clusterID);
            %                     if clusterSize <= maxClusterSize
            %                         numValidClusters = numValidClusters + 1;
            %                     end
            %                 end
            % 
            %                 %                 if numValidClusters < expectedClusters
            %                 %                     continue
            %                 %                 end
            %         
            %                 % Compute distances to template electrodes
            %                 distances = pdist2(selectedPoints,templateElectrode,'euclidean')
            %                 avgDistance = mean(distances);
            % 
            %                 % Compute curvature statistics
            %                 avgCurvature = mean(abs(selectedCurvatures));
            % 
            %                 % Normalize curvature and distances
            %                 normalized_curvature = avgCurvature / maxCurvature; % Assuming maxCurvature > 0
            %                 normalized_distance = 1 - (avgDistance / maxDistance); % Assuming maxDistance > 0
            % 
            %                 % Compute overall score
            %                 score = weight_curvature * normalized_curvature + weight_distance * normalized_distance;
            %                 
            %                 % Assign score
            %                 scores(k) = score;
            %             end
            
            %             % Find the parameter set with the highest score
            %             [maxScore, bestIdx] = max(scores);
            %             if isinf(maxScore) || maxScore == -Inf
            %                 disp(['No valid parameter combinations found for electrode ', num2str(i), '. Assigning template position.']);
            %                 all_nearest_points{i} = templateElectrode;
            %                 finalElectrodePositions(i, :) = templateElectrode;
            %                 continue;
            %             end
            % 
            %             bestThreshold = param_combinations(bestIdx, 1);
            %             bestEpsilon = param_combinations(bestIdx, 2);
            % 
            %             % Display the results
            %             disp(['Best Threshold: ', num2str(bestThreshold)]);
            %             disp(['Best Epsilon: ', num2str(bestEpsilon)]);
            % 
            %             % Use the best threshold and epsilon
            %             selectedVertices = (MC > bestThreshold); %| (MC < -bestThreshold);
            %             selectedPoints = scan_area.Points(selectedVertices, :);
            %             
            %             % Compute distances from the current template electrode to all segmented points
            %             if isempty(selectedPoints)
            %                 disp(['No valid points found after filtering for electrode ', num2str(i), '. Assigning template position.']);
            %                 all_nearest_points{i} = templateElectrode;
            %                 finalElectrodePositions(i, :) = templateElectrode;
            %                 continue
            %             end
            % 
            %             % Perform DBSCAN with best parameters
            %             [idx, ~] = dbscan(selectedPoints, bestEpsilon, minPoints);
            %             
            %             % Filter valid clusters
            %             uniqueClusters = unique(idx);
            %             validClusters = {};
            %             for clusterID = uniqueClusters'
            %                 if clusterID == -1 % Skip noise
            %                     continue;
            %                 end
            %                 % Get points in the current cluster
            %                 clusterPoints = selectedPoints(idx == clusterID, :);
            %                 if size(clusterPoints, 1) <= maxClusterSize
            %                     validClusters{end+1} = clusterPoints;
            %                 end
            %             end
            %             
            %             numValidClusters = numel(validClusters);
            %             disp(['Number of Valid Clusters: ', num2str(numValidClusters)]);
            %             
            %             % If no valid clusters found after filtering
            %             if numValidClusters == 0
            %                 disp(['No valid clusters found after filtering for electrode ', num2str(i), '. Assigning template position.']);
            %                 all_nearest_points{i} = templateElectrode;
            %                 finalElectrodePositions(i, :) = templateElectrode;
            %                 continue;
            %             end

            % Prepare colors
            %             colors = lines(numValidClusters);
            %             
            %             % Plot clusters
            %             figure(fig+2); clf;
            %             hold on;
            %             
            %             % Plot each valid cluster
            %             for i = 1:numValidClusters
            %                 clusterPoints = validClusters{i};
            %                 scatter3(clusterPoints(:, 1), clusterPoints(:, 2), clusterPoints(:, 3), ...
            %                          30, colors(i, :), 'filled');
            %             end
            %
            %             % Add labels and finalize
            %             xlabel('X');
            %             ylabel('Y');
            %             zlabel('Z');
            %             title(['Filtered Clusters with Threshold = ', num2str(bestThreshold), ...
            %                    ' and Epsilon = ', num2str(bestEpsilon)]);
            %             axis equal;
            %             grid on;
            %             hold off;
            %             disp(' ');
            
            % ==========================
            % Selection of the nearest point to electrode
            % ==========================
            %             segmentedPoints = selectedPoints;
            %             
            %             % Compute distances from the current template electrode to all segmented points
            %             if size(segmentedPoints,1) == 0 || size(templateElectrode,1) == 0
            %                 all_nearest_points{i} = templateElectrode;
            %                 finalElectrodePositions(i, :) = templateElectrode;
            %                 continue
            %             end
            %             dists = vecnorm(segmentedPoints - templateElectrode, 2, 2);
            %             
            %             % Find the nearest point
            %             [minDist, nearestIdx] = min(dists);
            %             % Store the nearest point, index, and distance
            %             nearestPoints = segmentedPoints(nearestIdx, :);
            %             all_nearest_points{i} = nearestPoints;
            %             nearestIndices = nearestIdx;
            %             distances = minDist;
            %             
            %             figure(fig+3); clf; 
            %             % Plot template electrodes
            %             scatter3(templateElectrode(:,1), templateElectrode(:, 2), templateElectrode(:, 3), ...
            %                      100, 'r', 'filled', 'DisplayName', 'Template Electrodes');
            %             hold on;
            %             % Plot segmented points
            %             scatter3(segmentedPoints(:, 1), segmentedPoints(:, 2), segmentedPoints(:, 3), ...
            %                      30, 'b', 'filled', 'DisplayName', 'Segmented Points');
            %             
            %             % Plot nearest points
            %             scatter3(nearestPoints(:, 1), nearestPoints(:, 2), nearestPoints(:, 3), ...
            %                      100, 'g', 'filled', 'DisplayName', 'Nearest Points');
            %             
            %             % Draw lines connecting template electrodes to their nearest points
            %             line([templateElectrode(:, 1), nearestPoints(:, 1)], ...
            %                  [templateElectrode(:, 2), nearestPoints(:, 2)], ...
            %                  [templateElectrode(:, 3), nearestPoints(:, 3)], 'Color', 'k');
            %             
            %             % Add labels and finalize
            %             xlabel('X');
            %             ylabel('Y');
            %             zlabel('Z');
            %             title('Template Electrodes and Nearest Points');
            %             axis equal;
            %             grid on;
            %             hold off;
            %
            %             figure(fig+4);
            %             pcshow(nearestPoints,'r','MarkerSize',600)
            %             hold on;

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
        Y = templateElectrodes; 
    else
        Y = templateElectrodes(found_idx,:);
    end
    
    % Init full set of options %%%%%%%%%%
    opt.method='nonrigid'; % use nonrigid registration
    
    opt.beta=4;            % the width of Gaussian kernel (smoothness)
    opt.lambda=6;          % regularization weight
    
    opt.viz=1;              % show every iteration
    opt.outliers = 1 - (sum(found_idx) / size(found_idx,1))       % noise weight
    opt.fgt=0;              % do not use FGT (default)
    opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
    opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)
    
    opt.max_it=150;         % max number of iterations
    opt.tol=1e-10;           % tolerance
    
    [Transform, C]=cpd_register(X,Y, opt);
    C
    
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

    % Parameters
    averageSpacing = 20; % Average Spacing of electrodes. 
                         % Also Maximum allowable distance from template 
                         % electrode to nearest point
    maxDistance = ceil(averageSpacing * 0.75);
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
        distanceToNearest = norm(nearestPoints(i, :) - alignedElectrodes(i, :));
        
        % If the nearest point is outside the diameter, replace it with the template position
        if distanceToNearest > maxDistance
            fprintf('Electrode %s: No nearest point within diameter. Using template position.\n', labels{i});
            finalElectrodePositions(i, :) = alignedElectrodes(i, :);
            missalignedIdx(i) = true; 
            continue
        end

        % If the template point has higher MC value, choose it
        dists = vecnorm(scan.Points - alignedElectrodes(i, :), 2, 2);
        closeIdx = dists <= maxDistance; % Points within allowable distance
        closePoints = scan.Points(closeIdx, :);
        closeMCValues = MC_cap(closeIdx); % Corresponding MC values

        if isempty(closePoints)
            % Skip if no points are within the allowable distance
            continue
        end
        
        % Find the MC value of the point in nearestPoints(i, :)
        nearestPointIdx = find(ismember(scan.Points, nearestPoints(i, :), 'rows'), 1);
        if isempty(nearestPointIdx)
            fprintf('Electrode %s: Nearest point not found in scan points.\n', labels{i});
            continue;
        end
        nearestPointMC = MC_cap(nearestPointIdx); % MC value for nearestPoints(i, :)

        % Find the best MC value around the template position
        localDists = vecnorm(closePoints - alignedElectrodes(i, :), 2, 2); % Distances from template
        localIdx = localDists <= localRadius; % Indices of points within the local radius

        if any(localIdx)
            % Get the point with the highest MC value in the local radius
            [bestMCValue, bestLocalIdx] = max(closeMCValues(localIdx));
            bestMCPoint = closePoints(localIdx, :);
            bestMCPoint = bestMCPoint(bestLocalIdx, :); % Best point based on MC value
    
            % Compare the best MC value around the template to the MC value of nearestPoints(i, :)
            if bestMCValue > nearestPointMC
                fprintf('Electrode %s: Template position has higher MC value. Using template position.\n', labels{i});
                finalElectrodePositions(i, :) = bestMCPoint; % Update to the point with the highest MC
                betterMCpointIdx(i) = true; % Mark as better point found
            end
        end
    end

    disp('Completed correction.'); disp(' ')

    % =====================
    % Save Aligned Electrodes
    % =====================
    pcElectrodes_aligned = pointCloud(finalElectrodePositions);
    save("refined_electrodes.mat","pcElectrodes_aligned")

    visualize_correction()

    % =====================
    % Display the final electrode positions
    % =====================
    visualize_final_electrode_positions()

    pcFinalElectrodePositions = pointCloud(finalElectrodePositions);
    disp('Finished electrode refinement.'); disp(' ')



    % =================
    % Align Template Again
    % =================
    fig = fig + 1;
    figure(fig);

    X = finalElectrodePositions;
    Y = templateElectrodes;
    if includes_nasion
        Y = Y(1:end-1,:);
    end
    
    % Init full set of options %%%%%%%%%%
    opt.method='affine'; % use nonrigid registration
    
    % Nonrigid params
    % opt.beta=6;            % the width of Gaussian kernel (smoothness)
    % opt.lambda=10;          % regularization weight
    
    % Affine params
    opt.scale = 1;
    opt.rot = 1;


    opt.viz=1;              % show every iteration
    opt.outliers = 0.2;       % noise weight
    opt.fgt=0;              % do not use FGT (default)
    opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
    opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)
    
    opt.max_it=150;         % max number of iterations
    opt.tol=1e-10;           % tolerance
    
    [Transform, C]=cpd_register(X,Y, opt);
    
    alignedElectrodes = cpd_transform(templateElectrodes,Transform);

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
    % REFINE WITH TEMPLATE MATCHING
    % ====================







    



    % =======================================
    % NESTED FUNCTIONS
    % =======================================

    function [GC, MC] = computeCurvatures(scan_area, do_plot)
        % Setup parallel computation
        numChunks = 4; % Adjust based on available cores
        chunkIndices = floor(linspace(1, size(scan_area.Points, 1) + 1, numChunks + 1)); % Chunk boundaries
        
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
        
        % Histogram
        %         numBins = 50;
        %
        % Create a new figure for the histogram
        %         figure('Name', 'Histogram of Mean Curvature Values', 'NumberTitle', 'off');
        %         
        %         % Plot the histogram
        %         histogram(MC, numBins, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'k', ...
        %             'BinLimits', [-0.5, 0.5], ...
        %             'Normalization', 'count');
        %         
        %         % Customize the plot
        %         xlabel('Mean Curvature');
        %         ylabel('Frequency');
        %         title('Histogram of Mean Curvature Values');
        %         grid on;
        

        if do_plot
            % Plot
            tri = scan_area.ConnectivityList;
            x = scan_area.Points(:,1);
            y = scan_area.Points(:,2);
            z = scan_area.Points(:,3);

            fig = fig + 1;
            img = figure(fig); clf;
            set(img, 'Position', [100 100 1200 600]);
        
            subplot(1,2,1:2)
            hold on
            axis equal
            pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',MC,'FaceColor','interp','EdgeColor','none');
            caxis([-0.36 , 0.36])
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
                 100, 'm', 'filled', 'DisplayName', 'Missaligned Electrodes');
    
    
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
        pcshow(finalElectrodePositions(betterMCpointIdx,:),'y','MarkerSize',600) 
        
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
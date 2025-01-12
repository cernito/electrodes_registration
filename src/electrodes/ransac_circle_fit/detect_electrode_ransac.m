    % ====================
    % Find electrode (circle) shapes using RANSAC on prealigned positions
    % ====================
    
    selectedVertices = (MC_cap > 0.02); %(MC < -bestThreshold);
    selectedPoints = scan.Points(selectedVertices, :);
    
    [idx, ~] = dbscan(selectedPoints, 20, 30);
    
    % Filter valid clusters
    uniqueClusters = unique(idx);
    validClusters = {};
    validVertices = {};

    for clusterID = uniqueClusters'
        if clusterID == -1 % Skip noise
            continue;
        end
        
        % Get points in the current cluster
        clusterMask = (idx == clusterID);
        clusterPoints = selectedPoints(clusterMask, :);
        clusterVertices = find(selectedVertices);
        clusterVertices = clusterVertices(clusterMask);

        % Discard clusters with large diameter size
        validClusters{end+1} = clusterPoints;
        validVertices{end+1} = clusterVertices;
    
        % if size(clusterPoints, 1) <= maxClusterSize
        %   validClusters{end+1} = clusterPoints;
        % end
    end
    
    numValidClusters = numel(validClusters);
    validPoints = vertcat(validClusters{:});
    validVertices = vertcat(validVertices{:});
    disp(['Number of Valid Clusters: ', num2str(numValidClusters)]);
    
    
    % Plot clusters
    figure(2); clf;
    hold on;
    
    % Plot selected points
    keepMask = false(size(scan.Points,1),1);
    keepMask(validVertices) = true;
    clusters = filterTriangulation(scan,keepMask);
    
    
    %trisurf(clusters.ConnectivityList,clusters.Points(:,1),clusters.Points(:,2),clusters.Points(:,3),'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');
    tri = clusters.ConnectivityList;
    x = clusters.Points(:,1);
    y = clusters.Points(:,2);
    z = clusters.Points(:,3);
    %pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',MC(selectedVertices),'FaceColor','interp','EdgeColor','none');
    pcshow(clusters.Points,'MarkerSize',300)
    %caxis([-0.36 , 0.36])
    %colormap jet
    %colorbar
    
    % Add labels and finalize
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    rotate3d on
    grid on;
    hold off;
    disp(' ');


    % Vysegmentuj body s vysokou MC hodnotou
    %     if did_smooth
    %         th = 0;
    %     else
    %         th = 0.05;
    %     end
    th = 0.02;
    segmentedVertices = (MC_cap > th); %| (MC_cap < -th);
    segmentedPoints = scan.Points(segmentedVertices, :);
    segmentedMC = MC_cap(segmentedVertices);

    [idx, ~] = dbscan(segmentedPoints, 20, 30);

    % Filter valid clusters
    uniqueClusters = unique(idx);
    validClusters = {};
    validVertices = {};

    for clusterID = uniqueClusters'
        if clusterID == -1 % skip noise
            continue
        end
        
         % Get points in the current cluster
        clusterMask = (idx == clusterID);
        clusterPoints = segmentedPoints(clusterMask, :);
        clusterVertices = find(segmentedVertices);
        clusterVertices = clusterVertices(clusterMask);

        % Get points in the current cluster
        validClusters{end+1} = clusterPoints;
        validVertices{end+1} = clusterVertices;
    end

    segmentedPoints = vertcat(validClusters{:});
    segmentedVertices = vertcat(validVertices{:});

    % Předpřiprav KDTree
    treeAll = KDTreeSearcher(scan.Points);
    treeSegmented = KDTreeSearcher(segmentedPoints);

    % Initialize electrodes
    numElectrodes = size(alignedElectrodes,1);
    otherElecMask = true(numElectrodes,1);

    figure(13);
    fig = 13;

    figure;
    for e_idx=1:numElectrodes

        % 1. Podívat se na okolí kolem elektrody - vybrat mesh v okolí.
        %   ? Poloměr okolí vybrat na základě vzdálenosti od ?
        %   ? nejbližší elektordy.                           ?
        %
        % 2. Vysegmentovat body s MC hodnotou nad daný threshold
        %   ? Použít pevný threshold nebo iterovat ?
        %
        % 3. V dané oblasti pomocí RANSAC najít kruh s podmínkou na poloměr 
        %    mezi 4-6 a cílem o kruh s největší průměrnou MC hodnotou.
        %
        % 4. Střed daného kruhu je nová pozice elektordy.
        %   ? pokud má větší mean absolute curvature ?
        %
        currentElectrode = alignedElectrodes(e_idx,:);
        otherElecMask(e_idx) = false;
        otherElectrodes = alignedElectrodes(otherElecMask,:);
        otherElecMask(e_idx) = true;

        % Najdi vzdálenost k nejbližší elektrodě
        [~, distance] = knnsearch(otherElectrodes,currentElectrode,'K',1);
        searchRadius = distance * 1.5;
        
        if isempty(searchRadius)
            finalElectrodePositions(e_idx, :) = currentElectrode;
            continue
        end

        % Vysegmentuj okolí elektrody
        closestIdx = knnsearch(scan.Points,currentElectrode,'K',1);
        [searchAreaIdx, ~] = rangesearch(treeAll, scan.Points(closestIdx,:), searchRadius);
        allSearchPoints = scan.Points(searchAreaIdx{:},:);
        searchMC = MC_cap(searchAreaIdx{:});

        [segmentedSearchIdx, ~] = rangesearch(treeSegmented, currentElectrode, searchRadius);
        searchPoints = segmentedPoints(segmentedSearchIdx{:,:},:);

        if size(searchPoints,1) < 15
            finalElectrodePositions(e_idx, :) = currentElectrode;
            continue
        end
        
        treeSearch = KDTreeSearcher(allSearchPoints);
        electrodeRadius = 5;
        [electrodeAreaIdx, ~] = rangesearch(treeSearch, currentElectrode, electrodeRadius);
        electrodeMC = searchMC(electrodeAreaIdx{:});


        [points2D, basis, meanPoint] = project_3D_to_plane(searchPoints);
        

        % Score
        N = size(points2D,1);

        distances = dist(points2D, currentElectrode(1), currentElectrode(2), 5);
        inlier_mask = abs(distances) < 1;
        inlier_count = sum(inlier_mask);

        meanMC = mean(abs(electrodeMC));
        stdMC = std(abs(electrodeMC));
    
        % Scoring metric
        bestScore = (inlier_count / N) * 0.5 + meanMC * 0.5 - stdMC * 0.1;


        [x0, y0, r, score] = fit_circle_ransac_MC(points2D, 6000, 1, bestScore, treeSearch, searchMC, basis, meanPoint);

        if r == 0
            finalElectrodePositions(e_idx, :) = currentElectrode;
            continue
        end

        if score < bestScore
            finalElectrodePositions(e_idx, :) = currentElectrode;
        end
        
        % Převeď 2D střed zp zpět do 3D
        % Center in 2D plane coordinates
        center2D = [x0; y0];
        center3D = meanPoint' + basis * center2D;
        
        fig = fig + 1;
        figure(fig)
        plot_circle(x0, y0, r, 'g', 'ransac');
        plot(points2D(:, 1), points2D(:, 2), '.r');
        axis equal;

        figure(13)
        plot_circle_3D(center3D, r, basis, searchPoints);
        
        % Radius zůstává stejný
        radius = r;

        % Odstraň vybrané body z hledného prostoru
        distancesToCenter = sqrt(sum((searchPoints - center3D'.^2),2));
        insideCircleMask = distancesToCenter <= (radius * 1.4);
        segmentedPoints = segmentedPoints(~ismember(segmentedPoints, searchPoints(insideCircleMask,:),'rows'),:);
        treeSegmented = KDTreeSearcher(segmentedPoints);
        segmentedMC = segmentedMC(~ismember(segmentedPoints, searchPoints(insideCircleMask, :), 'rows'));

        finalElectrodePositions(e_idx, :) = center3D;


    end
    pcFinalElectrodePositions = pointCloud(finalElectrodePositions);

    % =======================================
    % NESTED FUNCTIONS
    % ======================================

    function [points2D, basis, meanPoint] = project_3D_to_plane(points)
        meanPoint = mean(points, 1);
        centeredPoints = points - meanPoint;
        [~,~,V] = svd(centeredPoints, 0);
        basis = V(:,[1 2]);
        points2D = centeredPoints * basis;
    end

    function plot_circle_3D(center3D, r, basis, points)
        theta = linspace(0,2*pi);
        circle = center3D + r * basis * [cos(theta); sin(theta)];
        
        plot3(points(:,1),points(:,2),points(:,3),'.');
        hold on
        plot3(circle(1,:),circle(2,:),circle(3,:),'r','LineWidth', 2);
        axis equal
        scatter3(center3D(1,:),center3D(2,:),center3D(3,:),100,'red','filled')
    end
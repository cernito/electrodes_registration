scan = load_stl
og_scan = scan;
%%
scan = og_scan;
%%
pScan.vertices = scan.Points;
pScan.faces = scan.ConnectivityList;

% Number of faces (triangles) in the scan
numFaces = size(scan.ConnectivityList, 1);


% Set a dynamic reduction factor based on the number of vertices
if numFaces < 250000
    reduceFactor = 0.9;      % Less reduction for smaller scans
elseif numFaces < 450000
    reduceFactor = 0.6;     % Stronger reduction for detailed scans
else
    reduceFactor = 0.4;     % Stronger reduction for high-resolution scans
end

reduced_pScan = reducepatch(pScan,0.90);
new_scan.Points = reduced_pScan.vertices;
new_scan.ConnectivityList = reduced_pScan.faces;

scan = new_scan;

%%
%scan = og_scan;
FVscan.vertices = scan.Points;
FVscan.faces = scan.ConnectivityList;
            
smoothFVscan = smoothpatch(FVscan, 1, 2);
            
% Update the scan with the smoothed mesh
scan = triangulation(smoothFVscan.faces, smoothFVscan.vertices);

[GC, MC] = computeCurvatures(scan, true);
%[GC, MC] = curvatures(scan.Points(:, 1), scan.Points(:, 2), scan.Points(:, 3), scan.ConnectivityList);


%%
mean(GC,1,"omitnan")
selectedVertices = (GC > 0.05); %((MC < -0.08) & (MC > -0.02));
selectedPoints = scan.Points(selectedVertices, :);

[idx, ~] = dbscan(selectedPoints, 5, 5);
%idx = true(size(selectedPoints,1),1);

% Filter valid clusters
uniqueClusters = unique(idx);
validClusters = {};
validVertices = {};
maxClusterSize = 300;


for clusterID = uniqueClusters'
    if clusterID == -1 % Skip noise
        continue;
    end
    
    % Get points in the current cluster
    clusterMask = (idx == clusterID);
    clusterPoints = selectedPoints(clusterMask, :);
    clusterVertices = find(selectedVertices);
    clusterVertices = clusterVertices(clusterMask);

    validClusters{end+1} = clusterPoints;
    validVertices{end+1} = clusterVertices;

    %     if size(clusterPoints, 1) <= maxClusterSize
    %         validClusters{end+1} = clusterPoints;
    %         validVertices{end+1} = clusterVertices;
    %     end
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
ptClusters = pointCloud(clusters.Points);
ptClusters = pcdenoise(ptClusters);
pcshow(validPoints,'MarkerSize',300)
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


%%

% Example points
points = brushedData;
ptPoints = pointCloud(points);
ptPoints = pcdenoise(ptPoints);
points = ptPoints.Location;

%%

% Assuming you have a point cloud matrix
maxIter = 5000;
distanceThreshold = 1;

% Perform plane projection and circle fitting
[circleCenter, radius] = find_circle_ransac_svd(points, maxIter, distanceThreshold);

% Compute eigenvectors for visualization
covMatrix = cov(points);
[eigenVectors, ~] = eig(covMatrix);
[~, sortIdx] = sort(diag(eigenVectors), 'descend');
eigenVectors = eigenVectors(:, sortIdx);

% Visualize the result
meanPoint = mean(points, 1);

% ==================

% Identify points inside the circle
distancesToCenter = sqrt(sum((points - circleCenter').^2, 2));
insideCircleMask = distancesToCenter <= radius*2;

% Remove points inside the circle
points = points(~insideCircleMask, :);

x0 = circleCenter(1);
y0 = circleCenter(2);
r = radius;

hold on
plot_circle(x0, y0, r, 'g', 'ransac');
axis equal;

% ==================

%visualize_point_cloud_plane(points, meanPoint, eigenVectors);





%% 
% COMPUTATION OF GC ON SMALL AREA

radius = 30;
centroid = [-31.6, 36.5, 68.25];
scan_area = segment_scan_in_radius(scan, centroid, radius);
% Compute curvature values
[GC_scan_area, MC_scan_area] = computeCurvature(scan_area);
plot_MC(scan_area, MC_scan_area);
plot_GC(scan_area, MC_scan_area);

%%
selectedVertices = (MC_scan_area > 0.03);
selectedPoints = scan_area.Points(selectedVertices, :);

figure;
pcshow(selectedPoints,'MarkerSize',300)
hold on 
pcshow(scan_area.Points)





















%%
meanPoint = mean(points, 1);
covMatrix = cov(points);

% Compute eigenvectors and eigenvalues
[eigenVectors, eigenValues] = eig(covMatrix);

% Sort eigenvalues in descending order
[~, sortIdx] = sort(diag(eigenValues), 'descend');
eigenVectors = eigenVectors(:, sortIdx);

% The eigenvector corresponding to the smallest eigenvalue is the normal to the best-fitting plane
normalVector = eigenVectors(:, end)';

% Project points onto the plane
% Projection matrix: P = I - nn^T (where n is the normal vector)
projectionMatrix = eye(3) - (normalVector' * normalVector);

% Project each point onto the plane
projectedPoints = (projectionMatrix * (points - meanPoint)')' + meanPoint;

% Transform projected points to 2D plane coordinates
% Use the two principal eigenvectors as the basis
basis1 = eigenVectors(:, 1)';
basis2 = eigenVectors(:, 2)';

% Project points onto these basis vectors
points2D = zeros(size(projectedPoints, 1), 2);
for i = 1:size(projectedPoints, 1)
    % Compute displacement from mean point
    displacement = projectedPoints(i, :) - meanPoint;
    
    % Project displacement onto basis vectors
    points2D(i, 1) = dot(displacement, basis1);
    points2D(i, 2) = dot(displacement, basis2);
end

figure; clf;
subplot(1,2,1);
[d e f] = fit_circle_nhom(points2D);
[x0 y0 r] = quad_to_center(d,e,f);
plot_circle(x0, y0, r, 'b', 'nhom');

[d e f] = fit_circle_hom(points2D);
[x0 y0 r] = quad_to_center(d,e,f);
plot_circle(x0, y0, r, 'g', 'hom');
plot(points2D(:, 1), points2D(:, 2), '.r');
axis equal;

% Přidat MC hodnoty aby RANSAC vybral kruh s největší průměrnou MC hodnotou
subplot(1,2,2);
[x0 y0 r] = fit_circle_ransac_MC(points2D,5000,0.5);
plot_circle(x0, y0, r, 'g', 'ransac')
plot(points2D(:, 1), points2D(:, 2), '.r');
axis equal;

%%

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

function [GC, MC] = computeCurvature(TR)
    [GC, MC] = curvatures(TR.Points(:, 1), TR.Points(:, 2), TR.Points(:, 3), TR.ConnectivityList);
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
    clim([-0.1, 0.1]);
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

function clipped_curvature = clip_curvature(curvature, lower_threshold, upper_threshold)
    % CLIP_CURVATURE Clips curvature values to within [lower_threshold, upper_threshold]
    %
    % Inputs:
    %   curvature        - Vector of curvature values (MC or GC)
    %   lower_threshold  - Minimum allowable curvature value
    %   upper_threshold  - Maximum allowable curvature value
    %
    % Output:
    %   clipped_curvature - Vector of clipped curvature values
    
    clipped_curvature = curvature;
    
    % Apply lower threshold
    clipped_curvature(clipped_curvature < lower_threshold) = lower_threshold;
    
    % Apply upper threshold
    clipped_curvature(clipped_curvature > upper_threshold) = upper_threshold;
end


function [GC, MC] = computeCurvatures(scan_area, do_plot)
        % Setup parallel computation
        numChunks = 2; % Adjust based on available cores
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
            
            img = figure(20); clf;
            set(img, 'Position', [100 100 1200 600]); 
            
            subplot(1,2,1)
            hold on
            axis equal
            pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',GC,'FaceColor','interp','EdgeColor','none');
            caxis([-0.1 , 0.1])
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
            caxis([-0.1 , 0.1])
            colormap jet
            colorbar
            xlabel('x')
            ylabel('y')
            title('Estimated MC');
            hold off
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
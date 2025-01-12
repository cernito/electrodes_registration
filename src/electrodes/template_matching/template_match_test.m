
mainPath = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile('C:\ČVUT\Bakalarka\electrodes_registration','src')));

scan = load_stl;
og_scan = scan;
%%
pc_og_scan = pointCloud(og_scan.Points);


% ===========================
% SMOOTH AND REDUCE SCAN
% ===========================

scan = reduce_and_smooth(scan);
[~, MC_scan] = computeCurvature(scan);
plot_MC(scan, MC_scan);


%% ===========================
% SMOOTH AND REDUCE TEMPLATE
% ============================
electrode = load_stl;
electrode = reduce_and_smooth(electrode);
[~, MC_electrode] = computeCurvature(electrode);
plot_MC(electrode, MC_electrode);


%% =================================
% ASSIGN CURVATURE VALUES TO POINT CLOUD COLOR
% ==================================

rgb_scan = MC_to_rgb(MC_scan);
pcScan = pointCloud(scan.Points,"Color",rgb_scan);

rgb_electrode = MC_to_rgb(MC_electrode);
pcElectrode = pointCloud(electrode.Points,"Color",rgb_electrode);

figure
pcshow(pcElectrode,'MarkerSize',600);
hold on
pcshow(pcScan)


%% ======================
% Segment head scan around electrode
% =======================

redu_scan = reduce_tri(og_scan,0.99);
pcRedu_scan = pointCloud(redu_scan.Points);

[pcMovedElectrode, centroid] = move_electrode_to_nearest_point(pcElectrode, redu_scan);


% Segment part of the scan around the electrode
radius = 25;
scan_area = segment_scan_in_radius(redu_scan, centroid, radius);

% Compute the curvatrure
[~, MC_scan_area] = computeCurvature(scan_area);
plot_MC(scan_area, MC_scan_area);

rgb_scan_area = MC_to_rgb(MC_scan_area);
pcScan_area = pointCloud(scan_area.Points,"Color",rgb_scan_area);

% Visualization of the moved electrode and scan area
figure;
pcshow(pcMovedElectrode.Location,'r', 'MarkerSize', 600);
hold on;
pcshow(pcScan_area, 'MarkerSize', 200);
title('Electrode Moved to Nearest Point on Scan');
%%

figure
pcshow(pcMovedElectrode,"MarkerSize",600)
hold on
pcshow(pcScan_area,"MarkerSize",200)


%% =========================
% PREALIGN TEMPLATE TO SCAN_AREA PLANE
% ==========================

pcAlignedElectrode = match_scan_plane(pcScan_area, pcMovedElectrode);
%pcAlignedElectrode = pcMovedElectrode;


%% ===================== 
% Registration - Color ICP
% ======================

moving = pcAlignedElectrode;
fixed = pcScan_area;

moving.Normal = pcnormals(moving,6);
fixed.Normal = pcnormals(fixed,6);

[~,movingReg] = pcregistericp(moving,fixed,"Metric","planeToPlaneWithColor",...
    'MaxIterations', 300, 'Tolerance', [0.001 0.005]);

figure
%pcshowpair(movingReg,fixed,"ViewPlane","YX","ColorSource","Color")
pcshow(movingReg,"MarkerSize",600)
hold on
pcshow(fixed)


%%

figure
pcshow(og_scan.Points)
hold on
pcshow(movingReg.Location,'r',"MarkerSize",600)
















function Mby3_rgb = MC_to_rgb(MC)
    % Given mean curvature values

    % Step 1: Normalize curvature values to range [0, 1]
    min_val = -0.3;  % Expected minimum
    max_val = 0.3;   % Expected maximum
    normalized_MC = (MC - min_val) / (max_val - min_val);
    normalized_MC = max(0, min(1, normalized_MC)); % Clip values to [0, 1]
    
    % Step 2: Map normalized values to the jet colormap
    num_colors = 256;                % Number of colors in colormap
    jet_colormap = jet(num_colors);  % Generate jet colormap
    color_indices = round(normalized_MC * (num_colors - 1)) + 1;
    rgb_values = jet_colormap(color_indices, :);
    
    % Step 3: Combine into M-by-3 RGB matrix
    Mby3_rgb = rgb_values;
    
    figure;
    % Optional: Visualize the color values for the curvature
    scatter3(1:length(MC), zeros(size(MC)), MC, ...
             50, Mby3_rgb, 'filled');
    xlabel('Point Index'); ylabel('Y'); zlabel('Curvature');
    title('Mean Curvature Colored by Jet Colormap');
    colorbar;
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
    clim([-0.36 , 0.36]);
    colormap jet
    colorbar
    xlabel('x')
    ylabel('y')
    title('Estimated MC');
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


function TR = smooth_tri(TR)
    % Prepare the mesh structure for smoothing
    fvTR.vertices = TR.Points;
    fvTR.faces = TR.ConnectivityList;
    
    % Perform smoothing using the smoothpatch function
    smoothFVtri = smoothpatch(fvTR, 1, 1, 0.2);
    TR = triangulation(smoothFVtri.faces, smoothFVtri.vertices);
end


function TR = reduce_and_smooth(TR)
    TR = reduce_tri(TR,60000);
    TR = smooth_tri(TR);
end


function pcAlignedElectrode = match_scan_plane(pcScan_area, pcElectrode)
    %% =========================
    % PREALIGN TEMPLATE TO SCAN_AREA PLANE (WITH FLIP CORRECTION)
    % ==========================
   
     %% Step 1: Compute Scan Area Normal
    scan_points = pcScan_area.Location;
    [coeff_scan, ~, ~] = pca(scan_points);
    normal_scan = coeff_scan(:,3); 

    % Ensure the scan normal points outward from the center (0,0,0)
    scan_centroid = mean(scan_points, 1);
    if dot(normal_scan, scan_centroid) < 0
        normal_scan = -normal_scan;
    end

    %% Step 2: Compute Electrode Bump Direction
    electrode_points = pcElectrode.Location;
    electrode_centroid = mean(electrode_points, 1);
    [coeff_electrode, ~, ~] = pca(electrode_points);
    bump_direction = coeff_electrode(:,3);

    % If bump faces inward relative to the scan normal, flip it
    if dot(bump_direction, normal_scan) < 0
        bump_direction = -bump_direction;
    end

    %% Step 3: Compute Rotation to Align Bump Direction to Scan Normal
    v = cross(bump_direction, normal_scan); % rotation axis
    s = norm(v);
    c = dot(bump_direction, normal_scan);

    if s > 1e-6
        vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        R = eye(3) + vx + (vx * vx)*((1 - c)/s^2);
    else
        R = eye(3); % Already aligned
    end

    %% Step 4: Rotate Around Electrode Centroid
    translated_points = electrode_points - electrode_centroid;
    rotated_points = (R * translated_points')';
    aligned_points = rotated_points + electrode_centroid;

    % Step 7: Update the electrode point cloud
    pcAlignedElectrode = pointCloud(aligned_points, 'Color', pcElectrode.Color);
    
    % Visualization of the pre-aligned electrode and scan area
    figure;
    pcshow(pcAlignedElectrode.Location, 'r', 'MarkerSize', 600);
    hold on;
    pcshow(pcScan_area, 'MarkerSize', 200);
    title('Pre-Aligned Electrode and Scan Area');
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
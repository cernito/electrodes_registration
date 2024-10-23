% Ask the user to pick a filtering method
disp(' ');
disp('USER INPUT:');
disp('Choose a filtering method:');
disp('1: Filter points below a plane');
disp('2: Filter points outside a sphere');
method = input('Enter 1 or 2: ');

ptCloud = alignedToSphere_rotated;
figure(68)
pcshow(ptCloud)

% Downsample both point clouds for faster plotting (adjust gridStep to control the subsampling)
gridStep = 3;  % Adjust this value; larger values reduce the number of points
subsampledPtCloud = pcdownsample(ptCloud, 'gridAverage', gridStep);  % Subsample the moving point cloud

% Apply the corresponding filtering method
switch method
    case 1
        filteredCloud = filterAbovePlane(subsampledPtCloud, ptCloud);  % Call plane filtering function
    case 2
        filteredCloud = filterOutsideSphere(subsampledPtCloud, ptCloud);  % Call sphere filtering function
    otherwise
        disp('Invalid selection. Please enter 1 or 2.');
        return;
end



% --- Functions ---

function filteredCloud = filterAbovePlane(sampledPtCloud, originalPtCloud)
    % Function to allow the user to filter points above a Z-coordinate plane

    % Display the point cloud
    figure;
    pcshow(sampledPtCloud.Location, 'MarkerSize', 50);
    title('Select a Z-value to filter points below');
    axis equal;
    xlabel('X')
    ylabel('Y')
    zlabel('Z')

    % Ask the user to input the Z-coordinate value of the plane
    disp('1.1: Select a Z-value to filter points below');
    z_plane_below = input('Enter the Z-coordinate of the plane: ');

    % Filter out points above the Z-plane
    points1 = sampledPtCloud.Location;
    belowPlaneIdx1 = points1(:, 3) >= z_plane_below;  % Indices of points below or on the plane
    filteredPoints1 = points1(belowPlaneIdx1, :);

    % Display the filtered point cloud
    figure;
    pcshow(filteredPoints1, 'MarkerSize', 50);
    title('Select a Z-value to filter points above');
    axis equal;
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    % Ask the user to input the Z-coordinate value of the plane
    disp('1.2: Select a Z-value to filter points above');
    z_plane_above = input('Enter the Z-coordinate of the plane: ');

    % Filter out points below the Z-plane
    points1 = filteredPoints1;
    abovePlaneIdx1 = points1(:, 3) <= z_plane_above;  % Indices of points below or on the plane
    filteredPoints1 = points1(abovePlaneIdx1,:);
    
    % Display the filtered point cloud
    figure;
    pcshow(filteredPoints1, 'MarkerSize', 50);
    title('Filtered point cloud');
    axis equal;
    xlabel('X')
    ylabel('Y')
    zlabel('Z')

    % Filter the original pointCloud
    points = originalPtCloud.Location;
    belowPlaneIdx = points(:, 3) >= z_plane_below;  % Indices of points below or on the plane
    filteredPoints = points(belowPlaneIdx, :);
    % Filter out points below the Z-plane
    points = filteredPoints;
    abovePlaneIdx = points(:, 3) <= z_plane_above;  % Indices of points below or on the plane
    filteredPoints = points(abovePlaneIdx,:);
    filteredCloud = pointCloud(filteredPoints);
    figure(69)
    pcshow(filteredCloud);
end


function filteredCloud = filterOutsideSphere(sampledPtCloud, originalPtCloud)
    % Function to allow the user to filter points outside a defined sphere

    % Display the point cloud
    figure;
    pcshow(sampledPtCloud.Location, 'MarkerSize', 50);
    title('Select the center and radius of the sphere');
    axis equal;
    xlabel('X')
    ylabel('Y')
    zlabel('Z')

    % Ask the user to input the center and radius of the sphere
    center = input('Enter the center of the sphere as [x, y, z]: ');
    radius = input('Enter the radius of the sphere: ');

    % Filter out points outside the sphere
    points = sampledPtCloud.Location;
    distances = sqrt(sum((points - center).^2, 2));  % Compute distances from the center
    insideSphereIdx = distances <= radius;  % Indices of points inside the sphere
    filteredPoints = points(insideSphereIdx, :);

    % Display the filtered point cloud
    figure;
    pcshow(filteredPoints, 'MarkerSize', 50);
    title(['Points inside the sphere (radius: ', num2str(radius), ')']);
    axis equal;
    xlabel('X')
    ylabel('Y')
    zlabel('Z')

    % Filter the original pointCloud
    % Filter out points outside the sphere
    points = originalPtCloud.Location;
    distances = sqrt(sum((points - center).^2, 2));  % Compute distances from the center
    insideSphereIdx = distances <= radius;  % Indices of points inside the sphere
    filteredPoints = points(insideSphereIdx, :);

    % Create a new point cloud with filtered points
    filteredCloud = pointCloud(filteredPoints);

end

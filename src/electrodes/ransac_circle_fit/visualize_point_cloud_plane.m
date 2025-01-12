function visualize_point_cloud_plane(pointCloud, meanPoint, eigenVectors)
    % Visualize point cloud and its best-fitting plane
    % 
    % Inputs:
    %   pointCloud - Nx3 matrix of 3D points
    %   meanPoint - Mean point of the point cloud
    %   eigenVectors - Eigenvectors from PCA (3x3 matrix)

    % Create a new figure
    figure; clf;
    
    % Plot the original point cloud
    scatter3(pointCloud(:,1), pointCloud(:,2), pointCloud(:,3), 'b.');
    hold on;
    
    % Plot the mean point
    scatter3(meanPoint(1), meanPoint(2), meanPoint(3), 'ro', 'filled');
    
    % Create plane visualization
    % Use the first two eigenvectors to create a plane
    basis1 = eigenVectors(:,1)';
    basis2 = eigenVectors(:,2)';
    normalVector = eigenVectors(:,3)';
    
    % Compute the extent of the point cloud to scale the plane
    pointExtent = max(pointCloud) - min(pointCloud);
    maxExtent = max(pointExtent) * 1.2;  % Add 20% extra to fully cover the points
    
    % Create a grid of points in the plane
    [X,Y] = meshgrid(linspace(-maxExtent/2, maxExtent/2, 20));
    
    % Calculate plane points
    planePts = zeros(size(X,1), size(X,2), 3);
    
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            % Project plane points using basis vectors
            planePt = meanPoint + X(i,j)*basis1 + Y(i,j)*basis2;
            planePts(i,j,:) = planePt;
        end
    end
    
    % Plot the plane
    surf(planePts(:,:,1), planePts(:,:,2), planePts(:,:,3), ...
        'FaceColor', 'green', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Plot eigenvectors
    quiver3(meanPoint(1), meanPoint(2), meanPoint(3), ...
            basis1(1), basis1(2), basis1(3), maxExtent/2, 'r', 'LineWidth', 2);
    quiver3(meanPoint(1), meanPoint(2), meanPoint(3), ...
            basis2(1), basis2(2), basis2(3), maxExtent/2, 'g', 'LineWidth', 2);
    quiver3(meanPoint(1), meanPoint(2), meanPoint(3), ...
            normalVector(1), normalVector(2), normalVector(3), maxExtent/2, 'b', 'LineWidth', 2);
    
    title('Point Cloud with Best-Fitting Plane');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    legend('Points', 'Mean Point', 'Plane', 'Basis 1', 'Basis 2', 'Normal');
    
    axis equal;
    grid on;
    hold off;
end
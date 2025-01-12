function [center3D, radius, projectedPoints] = find_plane_circle_projection(pointCloud, maxIter, distanceThreshold)
    % Function to project 3D point cloud onto best-fitting plane and find circle
    % 
    % Inputs:
    %   pointCloud - Nx3 matrix of 3D points
    %   maxIter - Maximum iterations for RANSAC circle fitting
    %   distanceThreshold - Distance threshold for RANSAC
    %
    % Outputs:
    %   center3D - 3D coordinates of the circle center
    %   radius - Radius of the fitted circle
    %   projectedPoints - Points projected onto the plane

    % Perform Principal Component Analysis (PCA)
    % Compute mean and covariance of the point cloud
    meanPoint = mean(pointCloud, 1);
    covMatrix = cov(pointCloud);
    
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
    projectedPoints = (projectionMatrix * (pointCloud - meanPoint)')' + meanPoint;
    
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
    
    % Visualize the projections
    visualize_projected_points(points2D, pointCloud, projectedPoints);

    % Use provided RANSAC circle fitting function
    [x0, y0, r] = fit_circle_ransac(points2D, maxIter, distanceThreshold);
    
    % Reconstruct 3D circle center
    % Center in 2D plane coordinates
    center2D = [x0, y0];
    
    % Convert 2D center back to 3D
    center3D = meanPoint + ...
        center2D(1) * basis1 + ...
        center2D(2) * basis2;
    
    % Radius remains the same
    radius = r;
end

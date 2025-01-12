function [center3D, radius] = find_circle_ransac_svd(pointCloud, maxIter, distanceThreshold)
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

    % Project to plane SVD
    meanPoint = mean(pointCloud, 1);
    centeredPoints = pointCloud - meanPoint;
    [~,~,V] = svd(centeredPoints, 0);
    Q = V(:,[1 2]);
    points2D = centeredPoints * Q;

    % Use provided RANSAC circle fitting function
    [x0, y0, r] = fit_circle_ransac(points2D, maxIter, distanceThreshold);
    
    % Reconstruct 3D circle center
    % Center in 2D plane coordinates
    center2D = [x0; y0];

    center3D = meanPoint' + Q * center2D;
    theta = linspace(0,2*pi);
    circle = center3D + r * Q * [cos(theta); sin(theta)];
    
    % Radius remains the same
    radius = r;
    
    figure;
    plot3(pointCloud(:,1),pointCloud(:,2),pointCloud(:,3),'.');
    hold on
    plot3(circle(1,:),circle(2,:),circle(3,:),'r','LineWidth', 2);
    axis equal
    scatter3(center3D(1,:),center3D(2,:),center3D(3,:),100,'red','filled')


end

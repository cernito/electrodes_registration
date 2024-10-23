%%% Filter points outside the sphere
% Origin
origin = mean(pcSphere.Location);

% Maximum distance from the origin
maxDistance = 1.3*(z_limits(2) - z_limits(1))/2;  

% Calculate the Euclidean distances from the origin
points = alignedToSphere_filtered.Location;
distances = sqrt(sum((points - origin).^2, 2));

% Filter points where the distance is less than or equal to maxDistance
filteredPoints = points(distances <= maxDistance, :);
alignedToSphere_filtered = pointCloud(filteredPoints);

fig = fig + 1; 
figure(fig);
pcshowpair(pcSphere, alignedToSphere_filtered)
title('Filtered scan and sphere')

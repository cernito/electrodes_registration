
X = brushedData;

XC = mean(X,1);
Y=X-XC;
[~,~,V]=svd(Y,0);
Q = V(:,[1 2]); % basis of the plane
Y=Y*Q;

xc=Y(:,1);
yc=Y(:,2);
M=[xc.^2+yc.^2,-2*xc,-2*yc];
% Fit ellipse through (xc,yc)
P = M\ones(size(xc));

a=P(1);
P = P/a;
r=sqrt(P(2)^2+P(3)^2+1/a); % radius
xyzc = XC' + Q*P(2:3); % center

theta = linspace(0,2*pi);
c = xyzc + r*Q*[cos(theta); sin(theta)]; % fit circle

close all
plot3(X(:,1),X(:,2),X(:,3),'.');
hold on
plot3(c(1,:),c(2,:),c(3,:),'r','LineWidth', 2);
axis equal
scatter3(xyzc(1,:),xyzc(2,:),xyzc(3,:),100,'red','filled')


%
% meanPoint = mean(pointCloud, 1);
% centeredPoints = pointCloud - meanPoint;
% [~,~,V] = svd(centeredPoints, 0);
% Q = V(:,[1 2]);
% points2D = centeredPoints * Q;
% projectedPoints = points2D;
% 
% 
% 
% % Visualize the projections
% % visualize_projected_points(points2D, pointCloud, projectedPoints);
% 
% % Use provided RANSAC circle fitting function
% [x0, y0, r] = fit_circle_ransac(points2D, maxIter, distanceThreshold);
% 
% % Reconstruct 3D circle center
% % Center in 2D plane coordinates
% center2D = [x0, y0];
% 
% % Convert 2D center back to 3D
% %     center3D = meanPoint + ...
% %         center2D(1) * basis1 + ...
% %         center2D(2) * basis2;
% center3D = centeredPoints' + Q * [x0; y0];
%

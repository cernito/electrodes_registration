% Compute centroids of both models
centroid_scan = mean(pcHead_scan.Location);
centroid_mri = mean(pcMri_model.Location);

%%% Translate head_scan to match target mri centroid
%translation_vector = centroid_mri - centroid_scan;
%pcHead_scan = pcHead_scan + repmat(translation_vector, size(pcHead_scan,1),1);
pcHead_scan = pcHead_scan.Location - centroid_scan;
pcHead_scan = pointCloud(pcHead_scan);


[x,y,z]=sphere(50);

x_limits = pcMri_model.XLimits;
y_limits = pcMri_model.YLimits;
z_limits = pcMri_model.ZLimits;

%x = x*min(abs(pcMri_model.XLimits));
%y = y*min(abs(pcMri_model.YLimits));
%z = z*min(abs(pcMri_model.ZLimits));
x = x*(x_limits(2) - x_limits(1))/2;
y = y*(y_limits(2) - y_limits(1))/2;
z = 1.05*z*(z_limits(2) - z_limits(1))/2;

x_t = (x_limits(1)+x_limits(2))/2;
y_t = (y_limits(1)+y_limits(2))/2;
z_t = (z_limits(1)+z_limits(2))/2;

X = reshape(x+x_t,[],1);
Y = reshape(y+y_t,[],1);
Z = reshape(z+z_t,[],1);
pcSphere = pointCloud([X,Y,Z]);

fig = fig + 1; 
figure(fig);
pcshowpair(pcSphere, pcHead_scan)
title('Translated scan to centroid')

fig = fig + 1; 
figure(fig);
pcshowpair(pcSphere, pcMri_model)
title('Sphere and MRI')

%%% filtering
% minDistance = 2;
% [labels, numClusters] = pcsegdist(pcHead_scan, minDistance);
% 
% fig = fig + 1;
% figure(fig);
% pcshow(pcHead_scan.Location, labels)
% title('Segments from pcsegdist maxDist=3.5')
% 
% counts = histcounts(labels, numClusters);
% [sortedCounts, sortedIndices] = sort(counts, 'descend');
% % Check if there are at least two clusters
% if numClusters >= 2
%     largestClustersIdx = sortedIndices(1:2); % Indices of the two largest clusters
% else
%     error('Not enough clusters found');
% end
% 
% % Create a logical index for the two largest segments
% isInLargestClusters = (labels == largestClustersIdx(1)) | (labels == largestClustersIdx(2));
% 
% % Filter the point cloud to include only points from the two largest clusters
% pcHead_scan_filtered = select(pcHead_scan, isInLargestClusters);
% 
% fig = fig + 1; 
% figure(fig);
% pcshow(pcHead_scan_filtered)
% title('Filtered scan using pcsegdist maxDist=2')
pcHead_scan_filtered = pcHead_scan;

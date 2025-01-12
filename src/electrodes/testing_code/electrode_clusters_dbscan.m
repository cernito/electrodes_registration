%% Electrode detection

% full_file_path = get_user_file_path('*.stl', 'Select the scan of patients head');
% head_scan = stlread(full_file_path);
% all_points = head_scan.Points;
all_points = F1_area;

figure(998)
pcshow(all_points)

%pcloud = pcdownsample(pointCloud(all_points),"random",0.7);
%points = pcloud.Location;
%size(points)
points = downsample(all_points,2);

figure(999)
pcshow(points,'r')
hold on
pcshow(all_points)
%%
%points = downsample(all_points,2);

distance_matrix = squareform(pdist(points,'euclidean'));
%%

% Min-points: greater then dimension + 1 of input data
minPts = 10;

% Epsilon: select by k-distance graph
kD = pdist2(points,points,'euclidean','Smallest',minPts);

figure(100)
plot(sort(kD(end,:)));
title('k-distance graph')
xlabel(['Points sorted with ',num2str(minPts),'th nearest distances'])
ylabel([num2str(minPts),'th nearest distances'])
grid

%%

max_samples = 40;
epsilon = {15,10,5,3,1,0.1,0.01};
minimumPts = {3,5,10,15,20};
fig = 100;

num_electrodes = 128;
closest_to_num_electrodes = 128;
best = {0, 0, 0};

for j = 1:length(minimumPts)
    for i = 1:length(epsilon)
        eps = epsilon{i};
        minPts = minimumPts{j};
        [idx,corePts] = dbscan(distance_matrix, eps, minPts,'Distance','precomputed');
        
        % Extract cluster centroids
        [counts, uniqueClusters] = groupcounts(idx);
        length(uniqueClusters);
        validClusters = uniqueClusters((uniqueClusters ~= -1) & (counts <= max_samples));
        length(validClusters);
        
        cluster_centers = [];
        for k = 1:length(validClusters)
            cluster_points = points(idx == validClusters(k),:);
            cluster_centers = [cluster_centers; mean(cluster_points,1)];
        end
    
        fig = fig + 1;
        figure(fig); cla(fig);

        numGroups = length(validClusters);
        gscatter(points(:,1),points(:,2),idx,hsv(numGroups));
        title(['minPts=',num2str(minPts),' eps=',num2str(eps),' clusters=',num2str(length(validClusters))])

        if abs(num_electrodes - length(validClusters)) < closest_to_num_electrodes
            closest_to_num_electrodes = abs(num_electrodes - length(validClusters));
            best{1} = eps;
            best{2} = minPts;
            best{3} = length(validClusters);
        end
    end
end

disp(['best eps:',num2str(best{1})])
disp(['best minPts:',num2str(best{2})])
disp(['Num clusters:',num2str(best{3})])

% Extract cluster centroids
uniqueClusters = unique(idx);
electrodePositions = [];

for i = 1:length(uniqueClusters)
    if uniqueClusters(i) ~= -1 % Ignore noise points
        clusterPoints = points(idx == uniqueClusters(i), :);
        electrodePositions = [electrodePositions; mean(clusterPoints, 1)];
    end
end

electrodePC = pointCloud(electrodePositions);

figure(73)
pcshow(electrodePC);

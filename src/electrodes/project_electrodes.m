%% Euclidean projection of electodes onto the head
headPoints = pcMri_model.Location;
electrodePoints = pcElectrodes_aligned.Location;

% Create a KDTree for the head points
headKDTree = KDTreeSearcher(headPoints);

% Find nearest points on the head for each electrode point
[indices, distances] = knnsearch(headKDTree, electrodePoints);

% Project the electrode cap points onto the head points
projectedPoints = headPoints(indices, :);
projectedPoints = pointCloud(projectedPoints);

figure(72)
pcshow(pcMri_model.Location, [0 0.4470 0.7410])
hold on
pcshow(pcElectrodes_aligned);
pcshow(projectedPoints.Location, 'r', 'MarkerSize', 200);

labels = electrodes.labels;
labelOffset = 1;
labelPositions = projectedPoints.Location;

% Display electrode labels
electrodesPositions = pcElectrodes_aligned.Location;
for i = 1:size(electrodesPositions, 1)
    pos = electrodesPositions(i,:);
    offset = labelOffset * sign(pos);
        
    dr = 3;

    if pos(1) > 0
        offset(1) = offset(1) + dr;
    elseif pos(1) < 0
        offset(1) = offset(1) - dr;
    end

    if pos(2) > 0
        offset(2) = offset(2) + dr;
    elseif pos(2) < 0
        offset(2) = offset(2) - dr;
    end

    labelPositions(i,:) = pos + offset;
end

text(labelPositions(:,1), labelPositions(:,2), labelPositions(:,3), ...
    labels, "HorizontalAlignment", "center", "VerticalAlignment", "bottom", "Color", "cyan");

hold off

% ICP parameters
metric = "pointToPoint";
extrapolate = true;
inlierRatio = 1;
maxIterations = 400;
tolerance = [0.001 0.01];

% Performing ICP
[tform2,pcElectrodes_aligned] = pcregistericp(pcElectrodes, pcCap, ...
    Metric=metric, ...
    InlierRatio=inlierRatio, ...
    MaxIterations=maxIterations, ...
    Tolerance=tolerance, ...
    Extrapolate=extrapolate, ...
    Verbose=true);

figure(71); clf
pcshow(pcCap.Location, [.7 .7 .7])
hold on
pcshow(pcElectrodes_aligned.Location, 'r', 'MarkerSize', 200)

labels = electrodes.labels;
labelOffset = 1;
labelPositions = pcElectrodes_aligned.Location;

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
    labels, "HorizontalAlignment", "center", "VerticalAlignment", "bottom", "Color", "red");

hold on
pcshow(pcMri_model)
hold off

electrodes.pos = pcElectrodes_aligned.Location;
save('electrodes.mat', 'electrodes')

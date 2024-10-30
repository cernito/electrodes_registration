function pcHead_scan_filtered = prepare_scan_for_icp(pcMri_model, pcHead_scan)
    
    % Compute centroids of both models
    centroid_scan = mean(pcHead_scan.Location);
    centroid_mri = mean(pcMri_model.Location);
    
    %%% Translate head_scan to match target mri centroid
    %translation_vector = centroid_mri - centroid_scan;
    %pcHead_scan = pcHead_scan + repmat(translation_vector, size(pcHead_scan,1),1);
    pcHead_scan = pcHead_scan.Location - centroid_scan;
    pcHead_scan = pointCloud(pcHead_scan);
    
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

end

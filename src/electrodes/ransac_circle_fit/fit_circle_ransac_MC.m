function [x0, y0, r, bestScore] = fit_circle_ransac_MC(X, num_iter, threshold, bestScore, tree, searchMC, basis, meanPoint)
% function [x0 y0 r] = fit_circle_ransac(X, num_iter, threshold)
%
% INPUT: 
% X: n-by-2 vector
%    with data
% num_iter: number of RANSAC iterations
% threshold: maximal  distance of an inlier from the circumference
%
% OUTPUT: 
% cartesian coordinates of the circle

max_inliers = 0;
N = size(X,1);

x0 = 0;
y0 = 0;
r = 0;

for i=1:num_iter
    % Náhodně vyber 3 body
    sample_indices = randperm(N,3);
    sample_points = X(sample_indices,:);
    
    % Fitni kruh těmito body
    [d, e, f] = fit_circle_nhom(sample_points);
    [xi0, yi0, ri] = quad_to_center(d, e, f);
    
    % Vypočti počet inlierů
    distances = dist(X, xi0, yi0, ri);
    inlier_mask = abs(distances) < threshold;
    inlier_count = sum(inlier_mask);
    
    % Skip circles with unrealistic radii
    if ri < 4 || ri > 5
        continue;
    end

    center2D = [xi0; yi0];
    center3D = meanPoint' + basis * center2D;

    electrodeRadius = 5;
    [circleAreaIdx, ~] = rangesearch(tree, center3D', electrodeRadius);
    circleMC = searchMC(circleAreaIdx{:});

    % Compute mean and standard deviation of MC values
    meanMC = mean(abs(circleMC));
    stdMC = std(abs(circleMC));

    % Scoring metric
    score = (inlier_count / N) * 0.5 + meanMC * 0.5 - stdMC * 0.1;
    
    % Update the best circle based on mean MC value
    if score > bestScore
        bestScore = score;
        max_inliers = inlier_count;
        x0 = xi0;
        y0 = yi0;
        r = ri;
    end

end

end
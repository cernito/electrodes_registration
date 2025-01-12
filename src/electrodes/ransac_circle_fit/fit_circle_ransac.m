function [x0, y0, r] = fit_circle_ransac(X, num_iter, threshold)
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

for i=1:num_iter
    % Náhodně vyber 3 body
    sample_indices = randperm(N,3);
    sample_points = X(sample_indices,:);
   
    % Fitni kruh těmito body
    [d, e, f] = fit_circle_nhom(sample_points);   
    [xi0, yi0, ri] = quad_to_center(d, e, f);
    
    if ri <= 0.2 || ri >= 6.5
        continue
    end

    % Vypočti počet inlierů
    distances = dist(X, xi0, yi0, ri);
    inlier_count = sum(abs(distances) < threshold);
    
    % Zapamatuj si nejlepší parametry kružnice
    if inlier_count > max_inliers
        max_inliers = inlier_count;
        x0 = xi0;
        y0 = yi0;
        r = ri;
    end

end

end

function num_iter = calculate_ransac_iterations(p, w, m)
    % p: Desired confidence level (e.g., 0.99)
    % w: Probability of selecting an inlier
    % m: Minimum number of points to fit the model
    
    if w == 0
        error('The inlier probability (w) cannot be zero.');
    end

    num_iter = ceil(log(1 - p) / log(1 - w^m));
end
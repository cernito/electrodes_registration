function cdf = compute_cdf(data, intensity_levels)
    if ~exist('intensity_levels', 'var')
        intensity_levels = 256;
    end

    % Compute histogram
    histogram_counts = histcounts(data, linspace(0, 1, intensity_levels+1));

    % Normalize histogram to get probabilities
    pdf = histogram_counts / sum(histogram_counts);

    % Compute CDF
    cdf = cumsum(pdf);
end
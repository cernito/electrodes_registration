function MC_matched = match_hists(MC_scan, MC_template, intensity_levels)

    if ~exist('intensity_levels','var')
        intensity_levels = 256;
    end

    % Remove NaN values from both curvature arrays
    MC_scan_valid = MC_scan(~isnan(MC_scan));
    MC_template_valid = MC_template(~isnan(MC_template));

    % Normalize the input curvature values to [0, 1]
    MC_scan_norm = (MC_scan_valid - min(MC_scan_valid)) / (max(MC_scan_valid) - min(MC_scan_valid));
    MC_template_norm = (MC_template_valid - min(MC_template_valid)) / (max(MC_template_valid) - min(MC_template_valid));

    % Compute the CDFs for both scan and template curvatures
    cdf_scan = compute_cdf(MC_scan_norm, intensity_levels);
    cdf_template = compute_cdf(MC_template_norm, intensity_levels);

    % Create histogram matching lookup table
    matching_lut = zeros(intensity_levels, 1);
    for i=1:intensity_levels
        [~, intensity_val] = min(abs(cdf_scan(i) - cdf_template));
        matching_lut(i) = (intensity_val - 1) / (intensity_levels - 1);
    end

    % Match the histograms
    MC_scan_interp = interp1(linspace(0, 1, intensity_levels), matching_lut, MC_scan_norm, 'linear', 'extrap');
    
    % Rescale matched values back to the original MC_scan_range
    MC_matched_valid = MC_scan_interp * (max(MC_scan_valid) - min(MC_scan_valid)) + min(MC_scan_valid);
    
    MC_matched = NaN(size(MC_scan));
    MC_matched(~isnan(MC_scan)) = MC_matched_valid;

end
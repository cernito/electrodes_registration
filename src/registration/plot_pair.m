function plot_pair(pcSphere, pcHead_scan, pcMri_model)
    fig = 3; 
    figure(fig);
    pcshowpair(pcSphere, pcHead_scan)
    title('Translated scan to centroid')
    
    fig = fig + 1; 
    figure(fig);
    pcshowpair(pcSphere, pcMri_model)
    title('Sphere and MRI')
end


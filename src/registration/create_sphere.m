function pcSphere = create_sphere(pcMri_model)
    [x,y,z]=sphere(50);
    
    x_limits = pcMri_model.XLimits;
    y_limits = pcMri_model.YLimits;
    z_limits = pcMri_model.ZLimits;
    
    %x = x*min(abs(pcMri_model.XLimits));
    %y = y*min(abs(pcMri_model.YLimits));
    %z = z*min(abs(pcMri_model.ZLimits));
    x = x*(x_limits(2) - x_limits(1))/2;
    y = y*(y_limits(2) - y_limits(1))/2;
    z = 1.05*z*(z_limits(2) - z_limits(1))/2;
    
    x_t = (x_limits(1)+x_limits(2))/2;
    y_t = (y_limits(1)+y_limits(2))/2;
    z_t = (z_limits(1)+z_limits(2))/2;
    
    X = reshape(x+x_t,[],1);
    Y = reshape(y+y_t,[],1);
    Z = reshape(z+z_t,[],1);
    pcSphere = pointCloud([X,Y,Z]);
    
end


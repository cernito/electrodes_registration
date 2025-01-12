function [] = plot_circle(x0, y0, r, color, label)
	phi = 0:0.01:2*pi;
	x = r * cos(phi) + x0;
	y = r * sin(phi) + y0;
	plot(x, y, color, 'DisplayName', label, 'LineWidth', 2);
    hold on;


end

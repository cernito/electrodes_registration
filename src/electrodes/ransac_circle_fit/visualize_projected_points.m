function visualize_projected_points(points2D, originalPoints, projectedPoints)
    % Visualize projected points
    
    % Create a new figure with two subplots
    figure; clf;
    
    % First subplot - 2D Projection
    subplot(1,2,1);
    scatter(points2D(:,1), points2D(:,2), 'b.');
    title('2D Projected Points');
    xlabel('Basis 1 Projection');
    ylabel('Basis 2 Projection');
    axis equal;
    grid on;
    
    % Second subplot - 3D Original vs Projected Points
    subplot(1,2,2);
    % Plot original points
    scatter3(originalPoints(:,1), originalPoints(:,2), originalPoints(:,3), 'r.');
    hold on;
    % Plot projected points
    scatter3(projectedPoints(:,1), projectedPoints(:,2), projectedPoints(:,3), 'b.');
    
    title('Original vs Projected Points');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    legend('Original Points', 'Projected Points');
    axis equal;
    grid on;
    hold off;
end
function point_cloud = load_pointcloud()

    head_model = load_stl();
    point_cloud = pointCloud(head_model.Points);

end
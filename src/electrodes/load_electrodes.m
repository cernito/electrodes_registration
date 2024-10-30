file_path = get_user_file_path('*.elc', 'elc');
electrodes = elc_read(file_path);
labels = electrodes.labels;
Rz = [0 -1 0; 1 0 0; 0 0 1];
positions = electrodes.pos * Rz;
pcElectrodes = pointCloud(positions);

pcwrite(pcElectrodes, 'electrodes.pcd', 'Encoding', 'ascii');

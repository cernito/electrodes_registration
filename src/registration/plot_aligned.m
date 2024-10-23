%% PLOTTING
disp('Plotting final results...');

% scan and mri models in space
fig = fig + 1;
figure(fig); clf; trimesh(head_scan);
hold on
trimesh(mri_model)
title('Scan and mri models in space')

% translated scan and mri
fig = fig + 1;
figure(fig); clf
pcshowpair(pcHead_scan, pcMri_model)
title('Translated scan and mri')
%hold on
%pcshow(corresponding_target_cloud)
% hold on
% pcshow(pHead_scan.vertices)

% regisered scan on mri
fig = fig + 1;
figure(fig); clf
pcshowpair(ptCloudAligned, pcMri_model)
xlabel('x')
ylabel('y')
zlabel('z')
title('Registered scan')

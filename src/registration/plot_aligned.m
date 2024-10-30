%% PLOTTING
disp('Plotting final results...');

% scan and mri models in space
fig = 4;
figure(fig); clf; 
pcshow(pcHead_scan);
hold on
pcshow(pcMri_model);
title('Scan and mri models in space')

% translated scan and mri
fig = fig + 1;
figure(fig); clf
pcshowpair(pcHead_scan_translated, pcMri_model)
title('Translated scan and mri')
%hold on
%pcshow(corresponding_target_cloud)
% hold on
% pcshow(pHead_scan.vertices)

% regisered scan on mri
fig = fig + 1;
figure(fig); clf
pcshowpair(prealigned_scan, pcMri_model)
xlabel('x')
ylabel('y')
zlabel('z')
title('Registered scan')

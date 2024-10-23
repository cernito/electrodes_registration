figure(69); clf
grayColor = [.7 .7 .7];
pcshow(pcCap.Location, 'w')
hold on
pcshow(pcElectrodes.Location, 'r', 'MarkerSize', 200)
hold on 
pcshow(pcHead_scan)
hold off
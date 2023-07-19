clc; close all; clear all;
load('proj6_arterialFlowData_qfi.mat');

%mean magnitude image
mean_mag_m = mean(mag_3d, 3);
figure
imagesc(mean_mag_m)
colormap(gray)
axis image
axis off
title('Mean Magnitude of the image')

%creating body mask
bodyMask_m = zeros(256, 256);
max_pixel = max(mag_3d(:));
for row=1:256
    for col=1:256
        pixels = squeeze(mag_3d(row, col, :));
        if mean(pixels)>=0.1*max_pixel
            bodyMask_m(row, col) =1;
        end
    end
end
%displaying the body mask
figure
imagesc(bodyMask_m)
colormap(gray)
axis image
axis off
title("Head Mask")

vz_3d = venc*(phase_3d/pi*2);

% Showing a movie of velocity versus time. Red is flow toward head,
% blue is toward feet:
nTimes = length(time_v);
maxVz = max(vz_3d(:));
minVz = min(vz_3d(:));
figure
for timeIndex = 1:nTimes
    imagesc(bodyMask_m .* vz_3d(:, :, timeIndex))
    axis image
    axis off
    % Setting color limits to visualize slow and fast flow:
    set(gca, 'CLim', [minVz, maxVz]/2)
    drawnow
    m(timeIndex) = getframe;
end
nLoops = 2;
fps = 2; % Frames per second.
movie(m, nLoops, fps)

%identifying ascending aorta
maxVz_m = max(vz_3d, [], 3);
figure;
imagesc(maxVz_m .* bodyMask_m)
axis image
axis off
colorbar
title('Maximum Vz across time')
[aaMask_m, aaX_v, aaY_v] = roipoly;
line(aaX_v, aaY_v,'color', 'w', 'LineWidth', 1.5)

%identifying descending aorta
minVz_m = min(vz_3d, [], 3);
figure;
imagesc(minVz_m .* bodyMask_m)
axis image
axis off
colorbar
title('Minimum Vz across time')
[daMask_m, daX_v, daY_v] = roipoly;
line(daX_v, daY_v,'color', 'w', 'LineWidth', 1.5)

%calculating velocity
aaVp_m = zeros(256, 256);
daVp_m = zeros(256, 256);
nTimes = length(time_v);
aaVz_v = zeros(1, nTimes);
daVz_v = zeros(1, nTimes);
for timeIndex = 1:nTimes
    aa = aaMask_m.*squeeze(vz_3d(:,:,timeIndex));
    da = daMask_m.*squeeze(vz_3d(:,:,timeIndex));
    aaVz_v(timeIndex) = sum(aa(:))/sum(aaMask_m(:));
    daVz_v(timeIndex) = sum(da(:))/sum(daMask_m(:));
end

figure
plot(time_v/1000, aaVz_v, 'b:', time_v/1000, daVz_v, 'r-', 'LineWidth', 1.2);
title('Mean velocity in ascending and descending aorta over time')
xlabel('time [s]')
ylabel('Mean velocity [cm/s]')
box off
grid on

%measuring volume
aaVol = sum(aaMask_m(:))*(dx/10)*(dy/10)*mean(aaVz_v)*((time_v(end)-time_v(1))/1000);
daVol = sum(daMask_m(:))*dx/10*dy/10*mean(daVz_v)*(time_v(end)-time_v(1))/1000;

%calculating fraction of blood that goes from ascending aorta to descending
%aorta
fracDiff = 100 * (aaVol - abs(daVol)) / aaVol;

%to create velocity profile
all_dist = [];
[vmax, ind] = max(aaVz_v);
aaVp = aaMask_m.*squeeze(vz_3d(:,:,ind));
[~, z] = max(aaVp(:));
[xc, yc] = ind2sub([256 256],z);
for row = 1:256
    for col = 1:yc
        if aaMask_m(row, col)==1
            dist = sqrt((row-xc)^2 + (col-yc)^2);
            all_dist = [all_dist dist];
        end
    end
end
[B, I] = sort(all_dist);
Bnew1 = unique(B);
v = zeros(1, length(Bnew1));
for row = 1:256
    for col = 1:yc
        if aaMask_m(row, col)==1
            dist = sqrt((row-xc)^2 + (col-yc)^2);
            [Lia, Locb] = ismember(dist, Bnew1);
            v(1, Locb) = (v(1, Locb) + aaVp(row,col))/2;
        end
    end
end
figure
plot(Bnew1, v, 'bo');
hold on

for row = 1:256
    for col = yc:256
        if aaMask_m(row, col)==1
            dist = sqrt((row-xc)^2 + (col-yc)^2);
            all_dist = [all_dist dist];
        end
    end
end
[B, I] = sort(all_dist);
Bnew2 = unique(B);
v1 = zeros(1, length(Bnew2));
for row = 1:256
    for col = yc:256
        if aaMask_m(row, col)==1
            dist = sqrt((row-xc)^2 + (col-yc)^2);
            [Lia, Locb] = ismember(dist, Bnew2);
            v1(1, Locb) = (v1(1, Locb) + aaVp(row,col))/2;
        end
    end
end

plot(-flip(Bnew2), flip(v1), 'bo');
hold on

x= [-flip(Bnew2), Bnew1(1, 2:end)];
v_total = [flip(v1), v(1, 2:end)];
p = polyfit(x, v_total, 2);
y1 = polyval(p, x);
plot(x, y1, 'r', 'LineWidth', 2)
grid on



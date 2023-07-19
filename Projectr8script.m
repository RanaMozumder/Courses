%%
clc; close all; clear all;
load('proj8Data.mat'); %loading data

%%
time_v = cumsum(dTime_v); % Time at the completion of each image
nTimes = length(time_v);

%%
% Displaying PET (and MRI) images at each time point:
figure
bgImageMax = max(bgPetImage_3d(:));

% Making 3D arrays of identical pages to display MR images in grayscale:
bgMrColor_3d = cat(3, bgMrImage_m, bgMrImage_m, bgMrImage_m)/max(bgMrImage_m(:));
cbMrColor_3d = cat(3, cbMrImage_m, cbMrImage_m, cbMrImage_m)/max(cbMrImage_m(:));

for timeIndex = 1:nTimes
    % Showing PET image of basal ganglia:
    subplot(2, 2, 1)
    bgPetImage_m = squeeze(bgPetImage_3d(:, :, timeIndex));
    imagesc(bgPetImage_m)
    set(gca, 'CLim', [0, bgImageMax])
    axis image; axis off
    title(['Time index = ', num2str(timeIndex)])
    % Showing PET image of cerebellum:
    subplot(2, 2, 2)
    cbPetImage_m = squeeze(cbPetImage_3d(:, :, timeIndex));
    imagesc(cbPetImage_m)
    set(gca, 'CLim', [0, bgImageMax])
    axis image; axis off
    % Showing MR image of basal ganglia:
    subplot(2, 2, 3)
    image(bgMrColor_3d)
    axis image; axis off
    % Showing MR image of cerebellum:
    subplot(2, 2, 4)
    image(cbMrColor_3d)
    axis image; axis off
    pause(1)
end

%%
figure
imagesc(cbMrImage_m)
title('Define cerebellum')
[cbRoiMask_m, cb_x, cb_y] = roipoly;
line(cb_x, cb_y, 'Color', 'r', 'LineWidth', 3)

figure
imagesc(mean(cbPetImage_3d, 3))
title('Cerebellum in PET image')
hold on
line(cb_x, cb_y, 'Color', 'r', 'LineWidth', 3)

%%
for time=1:nTimes
    cbRoiMean_v (1, time) = sum(cbRoiMask_m.*squeeze(cbPetImage_3d(:,:,time)), 'all')...
        /sum(cbRoiMask_m, 'all');
end

%%
figure
plot(time_v, cbRoiMean_v, 'LineWidth', 3)
title('cbRoiMean_v versus time')
xlabel('Time')
ylabel('Mean Signal Intensity')

%%
figure
imagesc(mean(bgPetImage_3d, 3))
title('Define Basal Ganglia')
[bgRoiMask1_m, bg_x, bg_y] = roipoly;
line(bg_x, bg_y, 'Color', 'r', 'LineWidth', 3)
hold on
[bgRoiMask2_m, bg_x, bg_y] = roipoly;
line(bg_x, bg_y, 'Color', 'r', 'LineWidth', 3)
bgRoiMask_m = bgRoiMask1_m + bgRoiMask2_m;

%%
for time=1:nTimes
    bgRoiMean_v (1, time) = sum(bgRoiMask_m.*squeeze(bgPetImage_3d(:,:,time)), 'all')...
        /sum(bgRoiMask_m, 'all');
end

figure
plot(time_v, bgRoiMean_v, 'LineWidth', 3)
title('bgRoiMean_v versus time')
xlabel('Time')
ylabel('Mean Signal Intensity')

%%
for time=1:nTimes
    if time == 1
        bgRoiX_v(1, time) = 0;
        bgRoiY_v(1, time) = 0;
    else
        bgRoiX_v(1, time) = trapz(time_v(1, 1:time), cbRoiMean_v(1, 1:time))/bgRoiMean_v(1, time);
        bgRoiY_v(1, time) = trapz(time_v(1, 1:time), bgRoiMean_v(1, 1:time))/bgRoiMean_v(1, time);
    end
end

figure
plot(bgRoiX_v, bgRoiY_v, '-', bgRoiX_v, bgRoiY_v, '+')

title('Logan plot')

%%
nFitPoints = 4;
meanImage = mean(bgPetImage_3d, 3);
maxPix = max(meanImage(:));
mask = zeros(64, 64);
for row=1:64
    for col=1:64
        if mean(squeeze(bgPetImage_3d(row, col,:)))>=0.05*maxPix
            mask(row, col) = 1;
        end
    end
end
figure
imagesc(mask)


%%
dvr_m = zeros(64, 64);
for row=1:64
    for col=1:64
        if mask(row, col) == 1
            pixels = squeeze((bgPetImage_3d(row, col,:)))';
            for time=1:nTimes
                bgRoiX1_v(1, time) = trapz(time_v(1, 1:time), cbRoiMean_v(1, 1:time), 2)/pixels(1, time);
                bgRoiY1_v(1, time) = trapz(time_v(1, 1:time), pixels(1, 1:time), 2)/pixels(1, time);
            end
            X = bgRoiX1_v(1, (end-nFitPoints+1): end);
            Y = bgRoiY1_v(1, (end-nFitPoints+1): end);
            s = polyfit(X, Y, 1);
            % Check for legal values:
            if (s(1) > 0)
                dvr_m(row, col) = s(1); % DVR = slope.
            else
                dvr_m(row, col) = 0;
            end
        end
    end
end


%%
% Displaying DVR map:
figure
imagesc(dvr_m)
colorbar
axis image
axis off
title('Distribution Volume Ratio Map')
colorMap_m = colormap;
nColors = size(colorMap_m, 1);
% Showing thresholded DVR on MRI:
figure
% Calculating a mask where the DVR values
% exceed half their maximum value. Name the mask dvrMask_m:
dvrMask_m = zeros(64, 64);
max_pix = max(dvr_m(:));
for row=1:64
    for col = 1:64
        if dvr_m(row, col)>0.5*max_pix
            dvrMask_m(row, col) = 1;
        end
    end
end


% Rendering pixels inside dvrMask in color. Use same colormap as above:
colorIndex_v = 1 + round((nColors-1) * dvr_m(:) / max(dvr_m(:)));
redDvr_m = dvrMask_m .* reshape(colorMap_m(colorIndex_v, 1), size(dvrMask_m));
greenDvr_m = dvrMask_m .* reshape(colorMap_m(colorIndex_v, 2), size(dvrMask_m));
blueDvr_m = dvrMask_m .* reshape(colorMap_m(colorIndex_v, 3), size(dvrMask_m));
% Displaying the rest of the image in gray:
maskedMrImage_m = (1-dvrMask_m) .* bgMrImage_m / max(bgMrImage_m(:));
color_3d = cat(3, maskedMrImage_m + redDvr_m, maskedMrImage_m + greenDvr_m, maskedMrImage_m + blueDvr_m);
image(color_3d)
axis image
axis off
title('DVR map')
















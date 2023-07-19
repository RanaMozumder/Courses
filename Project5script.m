clc; close all; clear all;

%loading the data
load('proj5data_qfi.mat');

% % Making a movie:
% figure
% for index = 1:nTimes
%     imagesc(squeeze(image_3d(:, :, index)))
%     % Set the intensity scale based on the first image:
%     if (index == 1)
%         cLim_v = get(gca, 'CLim');
%     else
%         set(gca, 'CLim', cLim_v)
%     end
%     axis image
%     axis off
%     colormap(gray)
%     title(['Brain Images', num2str(index)])
%     drawnow
%     mov(index) = getframe;
% end
% fps = 4; % frames per second.
% nReps = 1; % number of repetitions.
% movie(mov, nReps, fps)

%creating a head mask
headMask_m = zeros(128, 128);
max_pixel = max(image_3d, [], 'all');
for row=1:128
    for col=1:128
        pixels = squeeze(image_3d(row, col, :));
        if mean(pixels)>=0.1*max_pixel
            headMask_m(row, col) =1;
        end
    end
end
% %displaying the head mask
% imagesc(headMask_m)
% colormap(gray)
% axis image
% axis off
% title("Head Mask")

%creating the time-of-minuimum map
timeOfMin_m = zeros(128, 128);
for row=1:128
    for col= 1:128
        if headMask_m(row, col)==1
            [minimum, I] = min(squeeze(image_3d(row, col, :)));
            timeOfMin_m(row,col) = I*tr;
        end
    end
end
%displaying the time-of-minimum map
imagesc(timeOfMin_m)
colorbar
colormap(gray)
axis image
axis off
title("Time of Minimum Map")

%selecting slow and control region
[slowRoiMask_m, x1_v, y1_v] = roipoly;
line(x1_v, y1_v,'color', 'y', 'LineWidth', 1)
hold on
[controlRoiMask_m, x1_v, y1_v] = roipoly;
line(x1_v, y1_v,'color', 'y', 'LineWidth', 1)
hold off
slowRoiMask_m = slowRoiMask_m .* headMask_m;
controlRoiMask_m = controlRoiMask_m .* headMask_m;

% % displaying the selected region masks
% figure
% imagesc(slowRoiMask_m + controlRoiMask_m)
% colorbar
% colormap(gray)
% axis image
% axis off
% title("Slow and Control ROI Mask")

%taking the mean signal value for each region after introducing the
%contrast agent
slowRoiMean_v = zeros(1, nTimes);
controlRoiMean_v = zeros(1, nTimes);
for timeIndex = 1:nTimes
    image_m = squeeze(image_3d(:, :, timeIndex));
    slowRoi_m = image_m .* slowRoiMask_m;
    controlRoi_m = image_m .* controlRoiMask_m;
    %Enter your own code here to calculate the mean signal intesity, S
    %in each ROI at the current time:
    slowRoiMean_v(timeIndex) = sum(slowRoi_m(:))/sum(slowRoiMask_m(:));
    controlRoiMean_v(timeIndex) = sum(controlRoi_m(:))/sum(controlRoiMask_m(:));
end

%displaying the mean signal intensity, S in selected regions
figure
time_v = tr * (0:(nTimes-1));
plot(time_v, slowRoiMean_v, 'r-', time_v, controlRoiMean_v, 'b:')
title('ROI mean signal')
ylabel('Signal Intensity')
xlabel('Time [s]')
legend('Slow flow region', 'Control region')

%measuring the baseline signal intensity, S0
baselineTime = input('Enter the duration of the baseline (in seconds): ');
baseIndex_v = find(time_v < baselineTime);
slowBaseSignal = mean(slowRoiMean_v(baseIndex_v));
controlBaseSignal = mean(controlRoiMean_v(baseIndex_v));

%calculating R2* change
slowR2_v = -log(slowRoiMean_v/slowBaseSignal)/te;
controlR2_v = -log(controlRoiMean_v/controlBaseSignal)/te;

%displaying R2* change curve w.r.t. time
figure
time_v = tr * (0:(nTimes-1));
plot(time_v, slowR2_v, 'r-', time_v, controlR2_v, 'b:')
title('R2* change in both ROIs')
ylabel('R2* change')
xlabel('Time [s]')
legend('Constrast agent conc. in Slow region', 'Constrast agent conc. in Control region')

%calculating relative cerebral blood volume
CBV_ratio = trapz(time_v, slowR2_v)/trapz(time_v, controlR2_v)








clc; close all;

load('proj3Data.mat');

%displaying all three images in the same figure
figure
subplot 131
imagesc(pdw_m);
colormap(gray)
axis off
axis image
title('Proton Density Weighted')

subplot 132
imagesc(t1w_m);
colormap(gray)
axis off
axis image
title('T1 Weighted')

subplot 133
imagesc(t2w_m);
colormap(gray)
axis off
axis image
title('T2 Weighted')

% %selecting training points
% nTrain = 10;
% figure
% imagesc(t2w_m)
% colormap(gray)
% axis image
% disp(['Click on’, num2str(nTrain), ‘ lesion points'])
% [x_v, y_v] = ginput(nTrain);
% row_v = round(y_v);
% col_v = round(x_v);
% 
% % %displaying the selected training points for lesion region
% % t2Max = max(t2w_m(:));
% % red_m = t2w_m/t2Max;
% % green_m = t2w_m/t2Max;
% % blue_m = t2w_m/t2Max;
% % for i=1:nTrain
% %     red_m(row_v(i), col_v(i)) = 1;
% %     green_m(row_v(i), col_v(i)) = 0;
% %     blue_m(row_v(i), col_v(i)) = 0;
% % end
% % color_3d = cat(3, red_m, green_m, blue_m);
% % figure
% % image(color_3d)
% % axis image
% % title('Point selected for training')
% 
% L_points = insertMarker(t2w_m/max(t2w_m(:)),[col_v row_v],'o','color','r', 'Size', 1);
% % figure
% % image(L_points)
% % axis image
% % axis off
% % title('Point selected for training')
% 
% %finding the intensity in PD and t2 weighted image for the selected points
% index_v = sub2ind(size(t2w_m), row_v, col_v);
% t2wl_v = t2w_m(index_v); % T2 weighted intensity for lesion points.
% pdwl_v = pdw_m(index_v); % PD weighted intensity for lesion points.
% 
% %mean intensity value for lesion region in PD and t2 weighted image
% t2wl = mean(t2wl_v);
% pdwl = mean(pdwl_v);
% 
% %displaying the selected training points for White Matter region
% nTrain = 10;
% figure
% imagesc(t2w_m)
% colormap(gray)
% axis image
% disp(['Click on’, num2str(nTrain), ‘ White Matter points'])
% [x_v, y_v] = ginput(nTrain);
% row_v = round(y_v);
% col_v = round(x_v);
% 
% WM_points = insertMarker(L_points,[col_v row_v],'+','color','g', 'Size', 1);
% % figure
% % image(WM_points)
% % axis image
% % axis off
% % title('Point selected for training')
% 
% index_v = sub2ind(size(t2w_m), row_v, col_v);
% t2ww_v = t2w_m(index_v); % T2 weighted intensity for WM points.
% pdww_v = pdw_m(index_v); % PD weighted intensity for Wm points.
% 
% t2ww = mean(t2ww_v);
% pdww = mean(pdww_v);
% 
% %displaying the selected training points for Gray Matter region
% nTrain = 10;
% figure
% imagesc(t2w_m)
% colormap(gray)
% axis image
% axis off
% disp(['Click on', num2str(nTrain), ' gray points'])
% [x_v, y_v] = ginput(nTrain);
% row_v = round(y_v);
% col_v = round(x_v);
% 
% GM_points = insertMarker(WM_points,[col_v row_v],...
%     'x','color','blue', 'Size', 1);
% % figure
% % image(GM_points)
% % axis image
% % title('Point selected for training')
% 
% index_v = sub2ind(size(t2w_m), row_v, col_v);
% t2wg_v = t2w_m(index_v); % T2 weighted intensity for gm points.
% pdwg_v = pdw_m(index_v); % PD weighted intensity for gm points.
% 
% t2wg = mean(t2wg_v);
% pdwg = mean(pdwg_v);
% 
% %displaying the selected training points for CSF region
% nTrain = 10;
% figure
% imagesc(t2w_m)
% colormap(gray)
% axis image
% disp(['Click on', num2str(nTrain), ' CSF points'])
% [x_v, y_v] = ginput(nTrain);
% row_v = round(y_v);
% col_v = round(x_v);
% 
% CSF_points = insertMarker(GM_points,[col_v row_v],...
%     'star','color','white', 'Size', 1);
% figure
% image(CSF_points)
% axis image
% axis off
% title('Point selected for training')
% 
% index_v = sub2ind(size(t2w_m), row_v, col_v);
% t2wc_v = t2w_m(index_v); % T2 weighted intensity for csf points.
% pdwc_v = pdw_m(index_v); % PD weighted intensity for csf points.
% 
% t2wc = mean(t2wc_v);
% pdwc = mean(pdwc_v);
% 
% figure
% plot(pdwl_v, t2wl_v, 'r.', pdww_v, t2ww_v, 'g+', ...
%     pdwg_v, t2wg_v, 'bx', pdwc_v, t2wc_v, 'k*')
% xlabel('PD-w intensity')
% ylabel('T2-w intensity')
% 
% pdMax = max(pdw_m(:));
% mask_m = (pdw_m > 0.1*pdMax);
% index2_v = find(mask_m(:));
% t2wHead_v = t2w_m(index2_v);
% pdwHead_v = pdw_m(index2_v);
% nHeadPixels = length(index2_v);
% 
% lesionMask_m = zeros(256, 256);
% wmMask_m = zeros(256, 256);
% gmMask_m = zeros(256, 256);
% csfMask_m = zeros(256, 256);
for pixel = 1:nHeadPixels
    %Find intensities for current pixel:
    t2w = t2wHead_v(pixel);
    pdw = pdwHead_v(pixel);
    %Insert code here to find the distance in feature space
    %between the current point at coordinates (pdw, t2w)
    %and the mean position of the training points for each
    %of the four pixel types (lesion, white matter, gray matter,
    %and cerebral spinal fluid).
    lesion_dist = sqrt((t2w-t2wl)^2 + (pdw-pdwl)^2);
    wm_dist = sqrt((t2w-t2ww)^2 + (pdw-pdww)^2);
    gm_dist = sqrt((t2w-t2wg)^2 + (pdw-pdwg)^2);
    csf_dist = sqrt((t2w-t2wc)^2 + (pdw-pdwc)^2);

    %Insert code here to find which pixel type the current
    %point is closest to.
    type = min([lesion_dist, wm_dist, gm_dist, csf_dist]);
    if type == lesion_dist
        lesionMask_m(ind2sub(size(t2w_m), index2_v(pixel))) = 1;
    elseif type == wm_dist
        wmMask_m(ind2sub(size(t2w_m), index2_v(pixel))) = 1;
    elseif type == gm_dist
        gmMask_m(ind2sub(size(t2w_m), index2_v(pixel))) = 1;
    else
        csfMask_m(ind2sub(size(t2w_m), index2_v(pixel))) = 1;
    end
    %Make the (row, col) element of the corresponding tissue
    %mask (lesionMask_m, etc) equal to 1. Note that the
    %current pixel%s position in the image is given by
    %[row, col] = ind2sub(size(t2w_m), index2_v(pixel));
end

figure
image(cat(3, lesionMask_m, wmMask_m, gmMask_m))
axis image
axis off

% disp('Define a polygon enclosing brain, excluding extracranial fat')
% skullMask_m = roipoly(pdw_m / pdMax);
% % Zero all pixels outside the skull:
% lesionMask_m = skullMask_m .* lesionMask_m;
% wmMask_m = skullMask_m .* wmMask_m;
% gmMask_m = skullMask_m .* gmMask_m;
% % Display improved segmentation map:
% figure
% image(cat(3, lesionMask_m, wmMask_m, gmMask_m))
% axis image
% axis off
% 
% n_lesion_pix = sum(lesionMask_m(:))
% n_wm_pix = sum(wmMask_m(:))
% n_gm_pix = sum(gmMask_m(:))
% n_csf_pix = sum(csfMask_m(:))
% 
% 
% 
% 
% 
% 

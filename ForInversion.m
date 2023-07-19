clc; close all;

load('proj1bData_QFI');
[nRows, nCols, nTe] = size(irImage_3d);

% plotting all the TI images in the same intensity scle
figure
subplot(3, 3, 1)
imagesc(squeeze(irImage_3d(:, :, 1)))
intLimits_v = get(gca, 'CLim');
axis image
axis off
colormap(gray)
title(['Inversion time = ', num2str(ti_v(1)), ' s'])

for index = 2:nTe
    subplot(3, 3, index)
    imagesc(squeeze(irImage_3d(:, :, index)), intLimits_v)
    axis image
    axis off
    colormap(gray)
    title(['Inversion time = ', num2str(ti_v(index)), ' s'])
end

% % creating a binary mask
% image_m = squeeze(irImage_3d(:, :, 7));
% mask_m = (image_m>0.2*max(image_m(:)));
% 
% figure
% imagesc(mask_m)
% colormap(gray)
% axis image
% axis off
% title('Binary Mask')
% 
% % creating T1 map
% t1_m = zeros(nRows, nCols);
% for row=1:nRows
%     for col=1:nCols
%         if (mask_m(row,col)==1)
%             Mz = squeeze(irImage_3d(row, col, :));
% 
%             numerator = 1- (Mz/m0_m(row,col));
%             ln = -real(log(numerator*0.5));
%             coeff_v = polyfit(ti_v, ln, 1);
%             
%             slope = coeff_v(1);
%         	t1 = 1 / slope;
%             t1_m(row, col) = t1;
%         end
%     end
% end
% 
% figure
% imagesc(t1_m)
% colormap(gray)
% colorbar
% 
% % creating residual map
% res=zeros(nRows,nCols);
% for row=1:nRows
%     for col = 1:nCols
%         signal_res_v = squeeze(irImage_3d(row, col, :));
%         res(row, col)=sqrt(norm(signal_res_v ...
%             - m0_m(row, col)*(1-2*exp(-ti_v./t1_m(row,col)))));
%     end
% end
% 
% figure
% imagesc(res)
% axis image
% axis off
% colormap(gray)
% title('Residual Map')
% 
% % showing fit for different regions
% figure
% subplot(1, 2, 1)
% t1Max = max(t1_m(:));
% red0_m = t1_m/t1Max;
% green0_m = t1_m/t1Max;
% blue0_m = t1_m/t1Max;
% color_3d = cat(3, red0_m, green0_m, blue0_m);
% image(color_3d)
% axis image
% title('T2 map: click on a point of interest')
% 
% nPoints = 50;
% contFlag = 1;
% while (contFlag == 1)
%     % Get position of mouse-click on image:
%     [x, y] = ginput(1);
%     row = round(y);
%     col = round(x);
%     % Exit loop if the mouse click is outside the image:
%     if (row < 1 || row > nRows || col < 1 || col > nCols)
%         contFlag = 0;
%         continue
%     end
%     red_m = red0_m;
%     green_m = green0_m;
%     blue_m = blue0_m;
%     % Show the position of the pixel in red:
%     red_m(row, col) = 1;
%     green_m(row, col) = 0;
%     blue_m(row, col) = 0;
%     color_3d = cat(3, red_m, green_m, blue_m);
%     image(color_3d)
%     axis image
%     title('T2 map: click on a point of interest')
%     % Show fit:
%     t1 = t1_m(row, col);
%     m0 = m0_m(row, col);
%     s_v = squeeze(irImage_3d(row, col, :));
%     
%     % Insert your code here to create an array of nPoints TE values from
%     % the minimum to maximum te:
%     tiFit_v = linspace(ti_v(1), ti_v(9), nPoints);
%     % Insert your code here to find the corresponding signal at each TE,
%     % using your estimates of T2 and S0:
%     sFit_v = m0*(1-2*exp(-tiFit_v./t1));
%     subplot(1, 2, 2)
%     plot(tiFit_v, sFit_v, ':', ti_v, s_v, '+')
%     title(['At row = ', num2str(row), ', col = ', num2str(col), ': T1 = ', ...
%         num2str(t1*1000), 'ms'])
% end
% 
% % creating mask for white matter
% figure
% image1_m = squeeze(irImage_3d(:,:,1));
% imagesc(image1_m);
% colormap(gray)
% axis off
% axis image
% [mask1_m, x1_v, y1_v] = roipoly;
% line(x1_v, y1_v,'color', 'y')
% 
% figure
% image1_m = squeeze(irImage_3d(:,:,1));
% imagesc(image1_m);
% colormap(gray)
% axis off
% axis image
% [mask2_m, x2_v, y2_v] = roipoly;
% line(x2_v, y2_v,'color', 'y')
% 
% figure
% image1_m = squeeze(irImage_3d(:,:,1));
% imagesc(image1_m);
% colormap(gray)
% axis off
% axis image
% [mask3_m, x3_v, y3_v] = roipoly;
% line(x3_v, y3_v,'color', 'y')
% 
% % combining all the masks
% mask_w_m=mask1_m + mask2_m +mask3_m;
% figure 
% imagesc(mask_w_m)
% axis image
% axis off
% colormap(gray)
% 
% % showing T1 map for white matter
% t1_w_m = zeros(nRows, nCols);
% for row=1:nRows
%     for col=1:nCols
%         if (mask_w_m(row,col)==1)
%             t1_w_m(row, col)=t1_m(row,col);
%         end
%     end
% end
% 
% figure 
% imagesc(t1_w_m)
% intLimits_v = get(gca, 'CLim');
% axis image
% axis off
% colormap(gray)
% 
% % histogram for gray matter
% figure
% histogram(t1_w_m)
% 
% % calculating mean and standard deviation for white matter
% w=1;
% for row=1:nRows
%     for col=1:nCols
%         if (mask_w_m(row,col)==1)
%             t1_w(w)=t1_m(row,col);
%             w=w+1;
%         end
%     end
% end
% 
% mean_white = mean(t1_w);
% StdDev_white = std(t1_w);
% 
% 
% % showing T1 map for gray matter
% figure
% image1_m = squeeze(irImage_3d(:,:,1));
% imagesc(image1_m);
% colormap(gray)
% axis off
% axis image
% [mask4_m, x4_v, y4_v] = roipoly;
% line(x4_v, y4_v,'color', 'y')
%  
% mask_g_m=mask4_m - mask1_m - mask2_m - mask3_m;
% figure 
% imagesc(mask_g_m)
% axis image
% axis off
% colormap(gray)
% 
% t1_g_m = zeros(nRows, nCols);
% for row=1:nRows
%     for col=1:nCols
%         if (mask_g_m(row,col)==1)
%             t1_g_m(row, col)=t1_m(row,col);
%         end
%     end
% end
% 
% figure 
% imagesc(t1_g_m)
% intLimits_v = get(gca, 'CLim');
% axis image
% axis off
% colormap(gray)
% histogram for gray matter
% figure
% histogram(t1_g_m)
% 
% calculating mean and standard deviation for gray matter
% g=1;
% for row=1:nRows
%     for col=1:nCols
%         if (mask_g_m(row,col)==1)
%             t1_g(g)=t1_m(row,col);
%             g=g+1;
%         end
%     end
% end
% 
% mean_gray = mean(t1_g);
% StdDev_gray = std(t1_g);
% 

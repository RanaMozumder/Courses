clc; close all;

load('proj1aData_QFI');
[nRows, nCols, nTe] = size(image_3d);

% plotting all the TI images in the same intensity scle
figure
subplot(2, 3, 1)
imagesc(squeeze(image_3d(:, :, 1)))
intLimits_v = get(gca, 'CLim');
axis image
axis off
colormap(gray)
title(['Echo time = ', num2str(te_v(1)*1000), ' ms'])

for index = 2:nTe
subplot(2, 3, index)
imagesc(squeeze(image_3d(:, :, index)), intLimits_v)
axis image
axis off
title(['Echo time = ', num2str(te_v(index)*1000), ' ms'])
end

% creating a binary mask
image_m = squeeze(image_3d(:, :, 1));
mask_m = (image_m>0.1*max(image_m(:)));

figure
imagesc(mask_m)
colormap(gray)
axis image
axis off
title('Binary Mask')

% creating T2 map
t2_m = zeros(nRows, nCols);
s0_m = zeros(nRows, nCols);
for row=1:nRows
    for col=1:nCols
        if (mask_m(row,col)==1)
            signal_v = squeeze(image_3d(row, col, :));
            coeff_v = polyfit(te_v, log(signal_v), 1);
            slope = coeff_v(1);
        	logS0 = coeff_v(2);	% Intercept.
        	t2 = -1 / slope;
            t2_m(row, col) = t2;
            s0_m(row, col) = exp(logS0);
        end
    end
end

for row=1:nRows
    for col = 1:nCols
        if t2_m(row, col)>0.2
            t2_m(row,col)=0.2;
        end
    end
end

figure
imagesc(t2_m, [0, 0.2])
colormap(gray)
colorbar

% creating residual map
res=zeros(nRows,nCols);
for row=1:nRows
    for col = 1:nCols
        %%%may need to change
        signal_res_v = squeeze(image_3d(row, col, :));
        res(row, col)=sqrt(norm(signal_res_v ...
            - s0_m(row,col)*exp(-te_v./t2_m(row,col))));
    end
end

figure
imagesc(res)
intLimits_v = get(gca, 'CLim');
axis image
axis off
colormap(gray)
title('Residual Map')

% showing fit for different regions
figure
subplot(1, 2, 1)
t2Max = max(t2_m(:));
red0_m = t2_m/t2Max;
green0_m = t2_m/t2Max;
blue0_m = t2_m/t2Max;
color_3d = cat(3, red0_m, green0_m, blue0_m);
image(color_3d)
axis image
title('T2 map: click on a point of interest')

nPoints = 50;
contFlag = 1;
while (contFlag == 1)
    % Get position of mouse-click on image:
    [x, y] = ginput(1);
    row = round(y);
    col = round(x);
    % Exit loop if the mouse click is outside the image:
    if (row < 1 || row > nRows || col < 1 || col > nCols)
        contFlag = 0;
        continue
    end
    red_m = red0_m;
    green_m = green0_m;
    blue_m = blue0_m;
    % Show the position of the pixel in red:
    red_m(row, col) = 1;
    green_m(row, col) = 0;
    blue_m(row, col) = 0;
    color_3d = cat(3, red_m, green_m, blue_m);
    image(color_3d)
    axis image
    title('T2 map: click on a point of interest')
    % Show fit:
    t2 = t2_m(row, col);
    s0 = s0_m(row, col);
    s_v = squeeze(image_3d(row, col, :));
    
    % Insert your code here to create an array of nPoints TE values from
    % the minimum to maximum te:
    teFit_v = linspace(te_v(1), te_v(5), nPoints);
    % Insert your code here to find the corresponding signal at each TE,
    % using your estimates of T2 and S0:
    sFit_v = s0*exp(-teFit_v/t2);
    subplot(1, 2, 2)
    plot(teFit_v, sFit_v, ':', te_v, s_v, '+')
    title(['At row = ', num2str(row), ', col = ', num2str(col), ': T2 = ', ...
        num2str(t2*1000), 'ms'])
end


% creating mask for white matter
figure
image1_m = squeeze(image_3d(:,:,1));
imagesc(image1_m);
colormap(gray)
axis off
axis image
[mask1_m, x1_v, y1_v] = roipoly;
line(x1_v, y1_v,'color', 'y')

figure
image1_m = squeeze(image_3d(:,:,1));
imagesc(image1_m);
colormap(gray)
axis off
axis image
[mask2_m, x2_v, y2_v] = roipoly;
line(x2_v, y2_v,'color', 'y')

mask_w_m = mask1_m - mask2_m;
figure 
imagesc(mask_w_m)
axis image
axis off
colormap(gray)

% showing T1 map for white matter
t2_w_m = zeros(nRows, nCols);
for row=1:nRows
    for col=1:nCols
        if (mask_w_m(row,col)==1)
            t2_w_m(row, col)=t2_m(row,col);
        end
    end
end

figure 
imagesc(t2_w_m)
intLimits_v = get(gca, 'CLim');
axis image
axis off
colormap(gray)

% histogram for gray matter
figure
histogram(t2_w_m)

% calculating mean and standard deviation for white matter
w=1;
for row=1:nRows
    for col=1:nCols
        if (mask_w_m(row,col)==1)
            t2_w(w)=t2_m(row,col);
            w=w+1;
        end
    end
end

mean_white = mean(t2_w);
StdDev_white = std(t2_w);

% creating mask for gray matter
figure
image1_m = squeeze(image_3d(:,:,1));
imagesc(image1_m);
colormap(gray)
axis off
axis image
[mask4_m, x4_v, y4_v] = roipoly;
line(x4_v, y4_v,'color', 'y')

mask_g_m = mask4_m - mask1_m;
figure 
imagesc(mask_g_m)
axis image
axis off
colormap(gray)

% showing T1 map for white matter
t2_g_m = zeros(nRows, nCols);
for row=1:nRows
    for col=1:nCols
        if (mask_g_m(row,col)==1)
            t2_g_m(row, col)=t2_m(row,col);
        end
    end
end

figure 
imagesc(t2_g_m, [0 0.2])
intLimits_v = get(gca, 'CLim');
axis image
axis off
colormap(gray)

% histogram for gray matter
figure
histogram(t2_g_m)

% calculating mean and standard deviation for white matter
w=1;
for row=1:nRows
    for col=1:nCols
        if (mask_g_m(row,col)==1)
            t2_w(w)=t2_m(row,col);
            w=w+1;
        end
    end
end

mean_gray = mean(t2_w);
StdDev_gray = std(t2_w);
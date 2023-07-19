load('proj4data.mat');
%displaying the original image
figure
imagesc(anat_m);
colormap(gray)

%creating a mask
mask_m = (anat_m>0.07*max(anat_m(:)));
figure
imagesc(mask_m)
colormap(gray)

%creating the FA matrix
fa_m = zeros(256, 256);
for row=1:256
    for col = 1:256
        if mask_m(row,col)==1
            lambda1 = eigValues_3d(row,col,1);
            lambda2 = eigValues_3d(row,col,2);
            lambda3 = eigValues_3d(row,col,3);
            lambda_bar = (lambda1 + lambda2 + lambda3)/3;
            numerator = (lambda1-lambda_bar)^2+(lambda2-lambda_bar)^2 ...
                +(lambda3-lambda_bar)^2;
            denominator = (lambda1)^2+(lambda2)^2+(lambda3)^2;
            fa_m(row, col) = sqrt(1.5*numerator/denominator);
        end
    end
end
%displaying the gray-scale FA map
figure
imagesc(fa_m)
colormap(gray)

% Display color-coded FA map:
red_m = fa_m .* abs(fastDiffVector_3d(:, :, 1));
green_m = fa_m .* abs(fastDiffVector_3d(:, :, 2));
blue_m = fa_m .* abs(fastDiffVector_3d(:, :, 3));
color_3d = cat(3, red_m, green_m, blue_m);
figure
imagesc(color_3d)
axis image
axis off
title('Red = R/L, Green = A/P, Blue = S/I')

%To define the seed point
figure
imagesc(anat_m)
axis image
colormap(gray)
hold on
disp('Define seed point...')
[x0, y0] = ginput(1);
hold on

x = x0;
y = y0;
x_v = x;
y_v = y;
stepFlag = 1;
cosAngle = 1;
stepSize = 1;
prevStep = [x_v,y_v];

while (stepFlag == 1)
    % To find the fast diffusion direction at the nearest
    % integer values of x and y:

    fast_v = squeeze(fastDiffVector_3d(round(y_v(end)),round(x_v(end)),1:2));

    % Break out of the while loop if the in-plane component of the
    % vector is too small:
    if (sum(fast_v.^2) < 0.5)
        disp('In-plane component of fast_v is too small')
        break
    end
    % To calculate cosAngle, the cosine of the angle between the previous step and
    % fast_v. If cosAngle is negative, reverse the direction of fast_v
    if length(x_v) > 1
        cosAngle = dot(prevStep,fast_v)/(norm(prevStep)*norm(fast_v));
        if cosAngle < 0
            fast_v = -fast_v;
        end
    end

    % Stepping a distance stepSize in the direction of fast_v
    % (i.e., update x_v and y_v):
    x = x+stepSize*fast_v(1);
    y = y+stepSize*fast_v(2);
    x_v = [x_v x];
    y_v = [y_v y];
    
    % line segment to the image to show the current step:
    line([x_v(end-1),x_v(end)],[y_v(end-1),y(end)], 'color', 'y')
    drawnow
    prevStep = fast_v;
    % setting stepFlag = 0 if abs(cosAngle) is too small:
    if abs(cosAngle)<0.5
        stepFlag = 0;
    end
    % End of while loop:
end
hold on
%OPPOSITE DIRECTION:

x2 = x0;
y2 = y0;
x2_v = x0;
y2_v = y0;
stepFlag = 1;
cosAngle2 = 1;
stepSize = 1;
prevStep2 = [x2_v,y2_v];
while(stepFlag == 1)

    fast2_v = -squeeze(fastDiffVector_3d(round(y2_v(end)),round(x2_v(end)),1:2));

    if (sum(fast2_v.^2) < 0.5)
        disp('In-plane component of fast_v is too small')
        break
    end

    if length(x2_v) > 1
        cosAngle2 = dot(prevStep2,fast2_v)/(norm(prevStep2)*norm(fast2_v));
        if cosAngle2 < 0
            fast2_v = -fast2_v;
        end
    end
    
    x2 = x2+stepSize*fast2_v(1);
    y2 = y2+stepSize*fast2_v(2);
    x2_v = [x2_v x2];
    y2_v = [y2_v y2];
    line([x2_v(end-1),x2_v(end)],[y2_v(end-1),y2(end)], 'color', 'y')
    drawnow
    prevStep2 = fast2_v;

    if abs(cosAngle2) < 0.5
        stepFlag = 0;
    end
end
hold off

%creating FA vector and matrix along the fiber path
fa = zeros(256, 256);
X = round([flip(x_v, 2), x2_v(1, 2:end)]);
Y = round([flip(y_v, 2), y2_v(1, 2:end)]);
for i=1:length(X)
    lambda1 = eigValues_3d(Y(1, i),X(1, i),1);
    lambda2 = eigValues_3d(Y(1, i),X(1, i),2);
    lambda3 = eigValues_3d(Y(1, i),X(1, i),3);
    lambda_bar = (lambda1 + lambda2 + lambda3)/3;
    numerator = (lambda1-lambda_bar)^2+(lambda2-lambda_bar)^2 ...
        +(lambda3-lambda_bar)^2;
    denominator = (lambda1)^2+(lambda2)^2+(lambda3)^2;
    fa_v(i) = sqrt(1.5*numerator/denominator);
    fa(Y(i),X(i))=fa_v(i);
end

%plotting FA vector vs position
figure
plot(1:length(X), fa_v)

%displaying the FA matrix
figure
imagesc(fa)
axis image
colormap(gray)
axis off
title("FA along fiber path")
colorbar










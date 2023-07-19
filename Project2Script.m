load('proj2Data')

%Make a movie:
figure
nFrames = 30;
for index = 1:nFrames
    imagesc(squeeze(image_3d(:, :, index)))
    %Setting the intensity scale based on the first image:
    if (index == 1)
        cLim_v = get(gca, 'CLim');
    else
        set(gca, 'CLim', cLim_v)
    end
    axis image
    axis off
    colormap(gray)
    title(['Cardiac phase: ', num2str(index)])
    drawnow
    mov(index) = getframe;
end
fps = 8; % frames per second.
nReps = 4; % number of repetitions.
movie(mov, nReps, fps)

% Defines the major and minor axes of the ellipse in each frame:
skip = 5; % Number of phases to skip.
nVols = length(1:skip:nFrames);
volume_v = zeros(1, nVols);
for index = 1:skip:nFrames
    imagesc(squeeze(image_3d(:, :, index)))
    set(gca, 'CLim', cLim_v)
    axis image
    axis off
    colormap(gray)
    hold on

    % Get axes:
    title(['Cardiac phase ', num2str(index), ' define major axis'])
    % Measures the major axis:
    [x_mj, y_mj] = ginput(2);
    row_mj = round(y_mj);
    col_mj = round(x_mj);
    d_mj = sqrt((row_mj(1)-row_mj(2))^2+(col_mj(1)-col_mj(2))^2);
    title(['Cardiac phase ', num2str(index), ' define minor axis'])
    % Measures the minor axis:
    [x_mn, y_mn] = ginput(2);
    row_mn = round(y_mn);
    col_mn = round(x_mn);
    d_mn = sqrt((row_mn(1)-row_mn(2))^2+(col_mn(1)-col_mn(2))^2);
    volIndex = (index-1)/skip + 1;
    % Calculates the LV volume:
    volume_v(volIndex) = (4*pi*d_mj*dx*d_mn*dx*d_mn*dx)/(3*8);
end

%calculates the ejection fraction
ejection_fraction = (max(volume_v)- min(volume_v))/max(volume_v);







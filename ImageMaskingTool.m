% Interactive script for masking bad regions in the X-ray diffraction    %%
% and scattering images produced at 14-ID-D beamline at APS              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all

%% Input

% Specify the path for image file you want to mask
image_path = 'D:\rimmerman_1702\insulin\protein\51_3ms_003.mccd';


% Specify the image center to help assess the quality of the mask
x_cen = 1990;   % X center of the image in pixels
y_cen = 1967;   % Y center of the image in pixels
detx = 365;     % sample to detector distance in mm
pixel_size = 34e4/3840*1e-3; % pixel size in mm
energy = 11.63;

% Do you want to create a new mask or use the previous one?
use_old_mask_flag = 1;     % 1 - use the old one, 0 - make a new one
old_mask_file = 'mask_insulin_24bm.h5';         % path for the old mask

% Do you want to save mask?
save_mask_flag = 0;         % 0 - do not save new mask, 1 - save new mask
save_mask_file = 'mask_insulin_24bm.h5';        % path where you want to create a new mask
% - you should not overwrite the old files
% - the file should have .h5 extension

%% Reading the image
image_data = imread(image_path);
image_data = double(image_data);

% Determine good scale for the image figure
[cts, bins] = hist(image_data(:),100);
zLim = [0 max(bins(cts>=0.1*max(cts)))];

%% Read or create mask

if use_old_mask_flag
    
    mask = readMask( old_mask_file );   
else 
    mask = zeros(size(image_data));

    % The edge pixels are masked by default. For most of detectors (including
    % Rayonix) edge pixels do not give adequate intensity and should not be
    % included into analysis.
    mask(1,:) = 1;
    mask(end,:) = 1;
    mask(:,1) = 1;
    mask(:,end) = 1;
end

%% Masking tool
% How to use this section:
% You have two knobs to control this section. The first is the Z-scale of the
% image, which help to see problematic regions. The second is the update_flag.
% If the flag is 1, then ths section will allow you to add a
% polygon mask to the image. If it is 0 then you will not be able to
% draw a polygon but only check if you are satisfied with the mask.
% The idea is as follows:
% 1) you run this section with "update_flag = 1"
% 2) draw a polygon to mask certain region. Double ckick the polygon so it is
% read in memory
% 3) close the figure and run this section again to draw the next polygon
% 4) repeat 1-3 until you are happy
% 5) put "update_flag = 0" to check the mask or just pass to the next
% section


% You can vary the image scale to better see which pixels have to be masked
Scale = 0.6;

% Update flag determines whether you want to draw another polygon or not.
% update_flag = 0 allows you to look at the image after all the masking
% has been performed 
update_flag = 0; % 0 - just shows the image; no polygon can be drawn; 1 - draws polygon

% plot the image
figure(1); clf;
imagesc(image_data,zLim(2)*[0.0 1.5]);
axis square
colormap('jet')

xlabel('X'); ylabel('Y')

% get the mask layer on top of the image in the black color
BLACK = cat(3, zeros(size(image_data)), zeros(size(image_data)), zeros(size(image_data)));
hold on
h = imagesc(BLACK);
hold off
set(h, 'AlphaData', mask);

if update_flag
    title({'Create polygon: ON','Double click the polygon once you have finished it'});
    mask_upd = roipoly;
    mask = mask | mask_upd;
else
    title('Create polygon: OFF');
end

%% Integrate the image with mask to see the pizel intensity distribution
q_grid_in = [];
distr_flag = 1;
% IntMax = zLim(2)*1.5;
IntMax = 4500;
clear azimuthalIntegrator
[ I, q, I_grid, I_distr ] = azimuthalIntegrator( image_data, mask, ...         % input image data
                                                 x_cen, y_cen, pixel_size, detx, energy, ... % input experimental parameters
                                                 q_grid_in, distr_flag, IntMax );                    % output space/grid

%% Check the mask quality
% This figure shows the distribution of the intensities per bin and the
% average curve on top of it. If you see some serious discontinuities in
% the behavior of these curves, redo the mask.

figure(2); clf;
hold on
imagesc(q,I_grid,I_distr')
colormap('gray');
colormap(1-colormap);

plot(q,I,'r-','linewidth',1)

box on;
xlabel('q, 1/A');
ylabel('intensity, cts')

xlim([q(1) q(end)])
ylim([I_grid(1) I_grid(end)])

%% Save the mask
if save_mask_flag
    h5create(save_mask_file,'/mask',size(single(mask)));
    h5write(save_mask_file,'/mask',single(mask));
end


function [ I, q_grid_out, I_grid, I_distr ] = azimuthalIntegrator( image_data, mask, ...    
                                            x_cen, y_cen, pixel_size, detx, energy, ... % input experimental parameters
                                            q_grid_in, distr_flag, IntMax )                        % output space/grid
%AZIMUTHALINTEGRATOR performs azimuthal integration of x-ray scattering
%patterns.
%In order to perform integration function requires a user to supply the
%following:
% image_data    - array containg 2D x-ray pattern
% mask          - array which masks undesired pixels
% x_cen, y_cen  - position of the beam center
% pixel_size    - size of the detector pixel in mm
% detx          - sample to detector distance in mm
% energy        - incident X-ray energy
% q_grid_in     - argument allowing to provide user defined grids
% of transfered momentum vector q. If user wants program to define the
% q-grid automatically the q_grid_in should be supplied as an empty array
% []. The proposed q-grid is returned as array contained in the q_grid_out.
% distr_flag    - allows user to obtain a distribution of pixel counts as a
% function of q. Typically, calculation of the distribution doubles the
% time required for integration of single image; therefore ...

persistent LUT Corrections q_grid

if isempty(LUT)
    [x_size, y_size] = size(image_data);
%     [X,Y] = meshgrid([1:x_size]+0.5-x_cen,[1:y_size]+0.5-y_cen);
    [X,Y] = meshgrid([1:x_size]-x_cen,[1:y_size]-y_cen);
    tth = atand(sqrt(X.^2+Y.^2)*pixel_size/detx);
    phi = atan2d(Y,X);
    q = 4*pi*energy/12.3987*sind(tth/2);
    P = 0.99;
    Corrections = (1/2*(1 + cosd(tth).^2 - P*cosd(2*phi).*sind(tth).^2)).* ...  % Polarization
                  (cosd(tth).^3);                                               % Geometry

    q = q(:);
    q = q(~mask(:));
    
    if isempty(q_grid_in)
        q_edges = linspace(min(q(:)), max(q(:))*(1+1e-6), 601);
        q_grid = q_edges(1:end-1) + diff(q_edges)/2;
    else
        q_grid = q_grid_in;
        q_edges = q_grid(1:end-1) + diff(q_grid)/2;
        q_edges = [q_grid(1)-abs(q_edges(1)-q_grid(1)); ...
                   q_edges; ...
                   q_grid(end)+abs(q_edges(end)-q_grid(end))];
    end
    
    LUT = defLUT(q, q_edges);
end
%% Integration of the image

image_data_corrected = image_data(:)./Corrections(:);
image_data_corrected = image_data_corrected(~mask(:));


I = zeros(length(q_grid),1);
if distr_flag
    I_edges = linspace(0, IntMax*(1+1e-6), round(IntMax));
    I_grid = I_edges(1:end-1) + diff(I_edges)/2;
    I_distr = zeros(length(q_grid), length(I_grid));
end

q_grid_out = q_grid;

integration_bench = tic;

for kk =1:length(q_grid)
    int_in_bin = image_data_corrected(LUT{kk});
    if ~isempty(int_in_bin)
%         Integration without dezingering:
%         I(kk) = mean( int_in_bin );
        
%         Integration with dezingering:
        int_mean = mean(int_in_bin);
        int_std  =  std(int_in_bin);
        int_incl = abs( int_in_bin - int_mean ) <= 3*int_std;
        I(kk) = mean( int_in_bin(int_incl) );
    end
    
    if distr_flag 
        dummy = histc(int_in_bin, I_edges);
        I_distr(kk,:) = dummy(1:end-1);
        I_distr(kk,:) = I_distr(kk,:)/max(I_distr(kk,:));
    end
    
end

disp(['Integration took ' num2str(toc(integration_bench)) ' s']);





    function LUT = defLUT(q, q_edges)
    [~, bin_idx] = histc(q, q_edges);
    LUT_bench = tic;

    bin_max = max(bin_idx);
    LUT = cell(size(q_edges));

    kk_0 = 0; kk_idx = 1;
    for kk =1:bin_max;
        dkk = (kk - kk_0)/bin_max*100;
        
        if dkk>=10;
            disp(['LUT creation progress: ' num2str(kk_idx*10) ' %']);
            kk_0 = kk; kk_idx = kk_idx + 1;
        end
        
        LUT{kk} = find(bin_idx==kk);
    end

    disp(['LUT creation took ' num2str(toc(LUT_bench)) ' s']);
    
    end

end


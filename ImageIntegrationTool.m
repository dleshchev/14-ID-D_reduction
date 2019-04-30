%%% Tool for integration of the images taken on 14-ID-D beamline at APS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
clear all
close all

%% INPUT
% paths for log files
log_paths = { ...
             'D:\rimmerman_1702\cytc\23.log', ...
             'D:\rimmerman_1702\cytc\24.log', ...
             'D:\rimmerman_1702\cytc\25.log', ...
             };
% paths for data. The size of this cell has to be the same as size of
% lig_paths cell.
data_paths = { ...
              'D:\rimmerman_1702\cytc\', ...
              'D:\rimmerman_1702\cytc\', ...
              'D:\rimmerman_1702\cytc\', ...
              };

% path for mask. Can be mccd, h5 and npy formats.
% mask_path = 'D:\work\Experiments\2016\12_2016_Insulin\Mask_rimmerman_316.npy';
mask_path = 'D:\matlab\CytC\mask_cytc_24bm.h5';
mask = readMask( mask_path );

energy = 11.63;     % Incident x-ray energy in keV
detx = 362.0;       % sample to detector distance in mm
x_cen = 1990.0;     % X beam position in pixels
y_cen = 1967.0;     % Y beam position in pixels
pixel_size = 340/3840; % pixel size in mm

% Desired q-grid for the curves
% q = 10.^linspace(log10(0.01),log10(4),600)';
q = linspace(0.005,4,1200)';

% q-range for normalization of the total scattering curves
%  qmin = 1.9; qmax = 2.1;
qmin = 1.475; qmax = 1.575;

% reference time delay for subtraction
t_ref = -5e-6; % in s

% Which reference OFF to subtract?
% ClosestRefFlag = 'Previous'; % previous before the current laser ON
% ClosestRefFlag = 'Closest'; % closest in lab time
ClosestRefFlag = 'MovingAverage'; NREF = 2; % Weighted average of the 
% previous and next laser OFF curves. NREF determines how many off curves 
% will be used. For example, NREF = 2 means that average of one previous 
% and one next laser OFF will be usefd. NREF MUST BE EVEN !!! 

% Saving options
save_data_flag = 1; % 1 - SAVES the file; 0 - does nothing
save_path = 'cytc_25C_18mgml_nts.h5'; % file for saving. has to have '.h5' extension

%% Reading, integrating and scaling the images
log_structure = [repmat('%s',1,4) repmat('%f',1,31)]; % auxiallry variable to aid with reading the log files
distr_flag = 0; % this is set to 0 because we don't need it in mass reduction of the data.

Data.total.s_raw = [];  % stores unscaled total scattering curves
Data.total.s     = [];  % stores scaled total scattering curves
Data.total.time  = [];  % stores time of creation of the image
Data.total.scan  = [];  % stores the name of the scan
Data.total.delay = [];  % stores time delay at which the curve was measured
          
n_scans = length(log_paths);
for ii = 1:n_scans
    
    [ log_data, nfiles ] = read_log_14idd( log_paths{ii} ); % get the log data
    
    % get the scan name string to log it along with other information
    logStr = strfind(log_paths{ii},'.log');
    slashes = strfind(log_paths{ii},'\');
    scan_str = log_paths{ii}(slashes(end)+1:logStr-1);
    
    % integrate stuff
    for kk = 1:nfiles
        tic;
        disp([scan_str ': ' num2str(round(ii/n_scans*100)) ' % / processed images in the scan: ' num2str(round(kk/nfiles*100)) ' %'])
        
        % awkward way of checking existance of the file
        im_path = [data_paths{ii} log_data.names{kk}];
        fid = fopen(im_path); fclose(fid);
        
        % if file does not exist then tell it
        if fid==-1
            disp(['file ' log_data.names{kk} ' is not found']);
        
        % if it exists then integrate it
        else
            im_data = double(imread(im_path));

            I = azimuthalIntegrator( im_data, mask, ...                          % input image data
                                     x_cen, y_cen, pixel_size, detx, energy, ... % input experimental parameters
                                     q, distr_flag, [] );                        % input q-grid
            I_scale = trapz( q(q>=qmin & q<=qmax), ...
                             I(q>=qmin & q<=qmax));

            % append integrated stuff and related attributes:
            Data.total.s_raw(:,end+1) = I;
            Data.total.s(:,end+1) = I/I_scale;
            Data.total.time(end+1) = log_data.times(kk);
            Data.total.scan{end+1} = scan_str;
            Data.total.delay(end+1) = log_data.delays(kk);
            
        end
        toc % this is just to tell you how much time you spend on each image
    end
end

%% calculate the differences

% get unique time delays and scan numbers. Note the logging the scan
% numbers works only with scans called as NUMBERS. In case of string/char
% scan names the scan names will not be saved
t = unique(Data.total.delay);
scannum = unique(Data.total.scan);

Data.diff.ds    = [];   % stores individual differential curves
Data.diff.delay = [];   % stores the time delay of the curves
Data.diff.scan  = [];   % stores the name of the scan
Data.diff.time  = [];   % stores the time when the LASER-ON image was collected

% the differential data is generated for each scan and delay separately and
% then stored together with appropriate flags.
for ii = 1:length(t)
    
    for kk = 1:length(scannum)
        
        % select the subset of the curves delay- and scan-wise
        subset_select = Data.total.delay == t(ii) & ...
                        strcmp(Data.total.scan, scannum(kk));
                    
        s_loc = Data.total.s(:,subset_select);
        time_loc = Data.total.time(subset_select);
              
        % if the delay is the reference one, then we will treat it slightly
        % differently. The way this statement works is a bit awkward
        % because of the double-precision of matlab.
        if abs(t(ii) - t_ref)<=1e-12
            
            ds_loc = s_loc(:,2:2:end) - s_loc(:,1:2:end-1);
            time_loc = time_loc(2:2:end);
        
        % otherwise we will take all the reference curvs along with the
        % laser-on curves. We look for the closest ones, but only in the
        % 'past' direction. This will make sure that each data collection
        % cycle is treated separately.
        else
        
            s_ref = Data.total.s(:, abs(Data.total.delay - t_ref) <=1e-12 & ...
                              strcmp(Data.total.scan, scannum(kk)));
        
            time_ref = Data.total.time( abs(Data.total.delay - t_ref) <=1e-12 & ...
                                     strcmp(Data.total.scan, scannum(kk)));
            % calculation of the time difference between laser-on and
            % reference curves:
            if strcmp(ClosestRefFlag, 'Previous')
                time_to_ref = bsxfun(@minus, time_loc, time_ref');
                % for the reference curves coming after the current laser-on
                % image we put large dummy numbers so they will be out of our
                % closest-partner search.
                time_to_ref(time_to_ref<0) = max(1e5, max(time_to_ref(time_to_ref>0)));
                [~, idx_ref] = min(time_to_ref,[],1);
                ds_loc = s_loc - s_ref(:,idx_ref);
            
            elseif strcmp(ClosestRefFlag, 'Closest')
                time_to_ref = abs(bsxfun(@minus, time_loc, time_ref'));
                [~, idx_ref] = min(time_to_ref,[],1);
                ds_loc = s_loc - s_ref(:,idx_ref);
            
            elseif strcmp(ClosestRefFlag, 'MovingAverage')
                ds_loc = zeros(length(q),length(time_loc));
                for jj = 1:length(time_loc)
                    
                    idx_pos = find(time_ref>time_loc(jj));
                    if ~(length(idx_pos)<fix(NREF/2))
                        idx_pos = idx_pos(1:fix(NREF/2));
                    end
                    
                    idx_neg = find(time_ref<time_loc(jj));
                    if ~(length(idx_neg)<fix(NREF/2))
                        idx_neg = idx_neg(length(idx_neg)-fix(NREF/2)+1:length(idx_neg));
                    end
                    
                    time_to_ref = abs(time_loc(jj)-time_ref([idx_neg idx_pos]));
                    ref_weight = 1./time_to_ref; ref_weight = ref_weight/sum(ref_weight);
                    ds_loc(:,jj) = s_loc(:,jj) - sum(bsxfun(@times, s_ref(:,[idx_neg idx_pos]), ref_weight),2);
                    
                end
            end
        end
        
        % append calcualted differentials
        Data.diff.ds =    [ Data.diff.ds ds_loc];
        Data.diff.delay = [ Data.diff.delay repmat(t(ii), size(ds_loc(1,:)))];
        Data.diff.scan  = [ Data.diff.scan  repmat(scannum(kk), size(ds_loc(1,:)))];
        Data.diff.time  = [ Data.diff.time  time_loc];

    end
end

%% plot the time delays.
% At this point it is important to go through all the time delays and look
% for obvious outliers. The following lines help to visualize the
% individual differences for selected time delay. The time delay selection
% is performed by changing 'nt' variable. You then identify the index of
% the bad curves and add them to idx_excl vector in the next section.

figure(1); clf;

% manually exclude some of the curves
idx_excl = []; % MovingAverage NREF = 2 case 


subplot(211)
hold on
nt = 17; % select the time delay

idx_subset = find(Data.diff.delay == t(nt));
curve_subset = Data.diff.ds(:,idx_subset);

plot(curve_subset);
% xlim([0.01, 3.3])
title(['t = ' num2str(t(nt))])
 
subplot(212)
plotyy(idx_subset,curve_subset(7,:),...
       idx_subset,curve_subset(400,:))
% plot(idx_subset,curve_subset(400,:))
grid on

% After you finish with throowing bad curves out, run this section one more
% time to make sure that variable 'idx_excl' is up to date and contains all
% the curves you want to throw out.

%% Averaging with outlier rejection
% here we will calculate the average of the curves. In addition we remove
% curves which are bad from Chavenet's criterion.
Data.diff.ds_av  = []; % stores averaged differential curves
Data.diff.ds_err = []; % stores error bars
Data.diff.t  = t;      % stores unique time delays

Data.diff.idx_incl = cell(length(t),1); % stores indexes of the curves used for average curve calculation

excl_mask = ones(size(Data.diff.ds(1,:)));
excl_mask(idx_excl) = 0;

% Criterion for the outlier rejection (Adapted Chauvenet's criterion).
% Values in a range of 0.03-0.15 usually work well. The small values
% correspond to 'hard' filtering (a lot of curves wil be rejected); the
% large values correspond to soft filtering.
Gamma = 0.25;

for ii = 1:length(t)
    % select a subset of curves:
    idx_subset = find(Data.diff.delay == t(ii) & excl_mask);
    curves_subset = Data.diff.ds(:,idx_subset);
    
    n_curves = length(curves_subset(1,:));
    
    % check each curve in the subset for 'outlier'-ness
    for kk = 1:n_curves
        diff_kk = curves_subset(:,kk);          % select a k-th curve
        diff_rest = curves_subset(:,:);         % select the rest of the subset
        diff_rest(:,kk) = [];                   % remove the k-th curve from the subset
        diff_rest_av = mean(diff_rest,2);       % get the average of the rest of the subset
        diff_rest_std = std(diff_rest,[],2);    % calculate the standard deviation of the dataset
        
        % figure of merit is the percentage of the points falling further
        % than 2 sigma away from the average curve.
        fom = nnz(abs(diff_kk-diff_rest_av)>2*diff_rest_std)/length(q);
        
        % Gamma determines the 'strictness' of the filtering.
        if fom<=Gamma;
            Data.diff.idx_incl{ii} = [Data.diff.idx_incl{ii} idx_subset(kk)];
        end
        
    end
    
    % After we have determined the indexes of the curves which have to be
    % included into the average calculation we can actually calculate the
    % average.
    Data.diff.ds_av(:,ii) = mean( Data.diff.ds(:,Data.diff.idx_incl{ii}), 2 );
    Data.diff.ds_err(:,ii) = std( Data.diff.ds(:, Data.diff.idx_incl{ii}), [], 2 )/sqrt(length( Data.diff.idx_incl{ii} ));
    disp(['[ t = ' num2str(t(ii)) ' s ] Outlier rejection results: ' num2str(length(Data.diff.idx_incl{ii})) ' out of ' num2str(n_curves) ' were accepted' ])
    
end

%% plot the results of the averaging along with individual curves

figure(2); clf; hold on;

ysh1 = 0.03;

for kk = 1:length(t)
    
    line([0,10],[0 0] - (kk-1)*ysh1,'color',[0.5 0.5 0.5])
%     plot(q, Data.diff.ds(:, Data.diff.delay == t(kk)) - (kk-1)*ysh1)
    plot(q, Data.diff.ds_av(:,kk) - (kk-1)*ysh1,'k.-')
    
end

% set(gca, 'xscale','log')
xlim([0.01 1])
box on;
xlabel('q, 1/A');
ylabel('\DeltaS(q,t)')


%% quick SVD to get an idea of what is going on
Q_sel = q>=0.02 & q<=0.6;
T_sel = 2:17;

DATA = Data.diff.ds_av(Q_sel, T_sel);

[uu,ss,vv] = svd(DATA);

figure(2); clf;

subplot(221)
plot(q(Q_sel),DATA)

subplot(222)
semilogy(diag(ss),'sq-')

n_cmp = 1:5;

subplot(223)
plot(q(Q_sel),uu(:,n_cmp))

subplot(224)
plot(Data.diff.t(T_sel),vv(:,n_cmp))
set(gca,'xscale','log')

%%

figure(3); clf;

plot(Data.diff.t(T_sel),vv(:,1),'sq-')
set(gca,'xscale','log')
xlim([5e-6 1e-1])
%% just look into individual average curves
figure(1); clf;
 
% semilogx(q,ds_av(:,3),'.-')
plot(q,Data.diff.ds_av(:,2),'.-')

xlim([0.001, 3.5])

%% time series with some integrals in different Q-regions

figure(5); clf;

subplot(121); hold on
ysh2 = 2e-2;

Color = colormap(cool(length(t)));
for ii = 1:length(t)
    
    plot([0.02 1],[0 0]-(ii-1)*ysh2,'color',[0.5 0.5 0.5])
    plot(q, Data.diff.ds_av(:,ii) - (ii-1)*ysh2,'.-','color',Color(ii,:));
    text(0.35,0.007-(ii-1)*ysh2,[num2str(Data.diff.t(ii)*1e3) ' ms'])
end

xlim([0.02 1])
set(gca,'xscale','log');
box on
xlabel('q, A^{-1}'); ylabel('\DeltaS(q,t)'); set(gca,'yticklabel',[])
title('Time series')

Qroi1 = [0.02 0.08]; S1 = sum(Data.diff.ds_av(q>=Qroi1(1) & q<=Qroi1(2),:),1);
Qroi2 = [0.08 0.15]; S2 = sum(Data.diff.ds_av(q>=Qroi2(1) & q<=Qroi2(2),:),1);
Qroi3 = [0.10 0.40]; S3 = sum(Data.diff.ds_av(q>=Qroi3(1) & q<=Qroi3(2),:),1);

subplot(222); hold on; h = [];
line([1e-12 1],[0 0],'color',[0.5 0.5 0.5])
h(1) = plot(Data.diff.t(2:end), -S1(2:end),'sq-');
% h(2) = plot(Data.diff.t(2:end),2*S2(2:end),'o-');
h(2) = plot(Data.diff.t(2:end),1.3*S3(2:end),'^-');

set(gca,'xscale','log');
xlabel('time, s'); ylabel('abs. int. signal, a.u.');
xlim([5e-7 200e-3]);
box on
% legend(h,[num2str(Qroi1(1)) '-' num2str(Qroi1(2)) ' A^{-1}'],...
...%          [num2str(Qroi2(1)) '-' num2str(Qroi2(2)) ' A^{-1}'],...
%          [num2str(Qroi3(1)) '-' num2str(Qroi3(2)) ' A^{-1}']);
legend(h,[num2str(Qroi1(1)) '-' num2str(Qroi1(2)) ' A^{-1}'],...
         [num2str(Qroi3(1)) '-' num2str(Qroi3(2)) ' A^{-1}']);
legend('boxoff')
title('Signal integration')

%% save stuff

if save_data_flag
    
    h5create(save_path,'/total/s_raw',size(Data.total.s_raw));
    h5write(save_path,'/total/s_raw',Data.total.s_raw);

    h5create(save_path,'/total/s',size(Data.total.s));
    h5write(save_path,'/total/s',Data.total.s);
    
    h5create(save_path,'/total/time',size(Data.total.time));
    h5write(save_path,'/total/time',Data.total.time);
    
    if ~isempty( str2num(cell2mat(Data.total.scan')) );
        h5create(save_path,'/total/scan', size( str2num(cell2mat(Data.total.scan')) ) );
        h5write(save_path,'/total/scan', str2num(cell2mat(Data.total.scan')) );
    end
    
    h5create(save_path,'/total/delay', size( Data.total.delay ) );
    h5write(save_path,'/total/delay',Data.total.delay);
    
    h5create(save_path,'/diff/ds', size( Data.diff.ds ) );
    h5write(save_path,'/diff/ds',Data.diff.ds);
    
    h5create(save_path,'/diff/time', size( Data.diff.time ) );
    h5write(save_path,'/diff/time', Data.diff.time);
    
    if ~isempty( str2num(cell2mat(Data.total.scan')) );
        h5create(save_path,'/diff/scan', size( str2num(cell2mat(Data.diff.scan')) ) );
        h5write(save_path,'/diff/scan', str2num(cell2mat(Data.diff.scan')) )
    end
    
    h5create(save_path,'/diff/delay', size( Data.diff.delay ) );
    h5write(save_path,'/diff/delay', Data.diff.delay);
    
    h5create(save_path,'/diff/ds_av', size( Data.diff.ds_av ) );
    h5write(save_path,'/diff/ds_av', Data.diff.ds_av );
    
    h5create(save_path,'/diff/ds_err', size( Data.diff.ds_err ) );
    h5write(save_path,'/diff/ds_err', Data.diff.ds_err);

    h5create(save_path,'/diff/t', size( Data.diff.t ) );
    h5write(save_path,'/diff/t', Data.diff.t);
    
    lens = cellfun('length',Data.diff.idx_incl);
    idx_incl = zeros(numel(lens),max(lens));
    for ii =1:length(lens)
        idx_incl(ii,[1:max(lens)]<=lens(ii)) = vertcat(Data.diff.idx_incl{ii});
    end
    
    h5create(save_path,'/diff/idx_incl', size( idx_incl ) );
    h5write(save_path,'/diff/idx_incl', idx_incl);
    
    h5create(save_path,'/q', size( q ) );
    h5write(save_path,'/q', q);
    
end



function [ log_data, nfiles ] = read_log_14idd( log_path )
%READ_LOG_14IDD Summary of this function goes here
%   Detailed explanation goes here

log_structure = [repmat('%s',1,4) repmat('%f',1,31)];

fid = fopen(log_path);
log_data = textscan(fid,log_structure,'commentstyle','#');
fclose(fid);

file_dates = log_data{1};
file_times = log_data{2};
file_names = log_data{3};
delays = log_data{4};

nfiles = length(file_names);

delays_num = zeros(size(delays));
file_times_s = zeros(size(delays));
for kk =1:nfiles
%     disp(delays{kk})
    delays_num(kk) = time_str2num(delays{kk});
    file_times_s(kk) = datenum([file_dates{kk} ' ' file_times{kk}],'dd-mmm-yy HH:MM:SS')*86400;
end

log_data = struct('names', {file_names}, ...
                  'times', file_times_s, ...
                  'delays', delays_num);


end


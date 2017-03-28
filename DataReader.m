function [  q, t, ds_av, ds_err, Data] =  DataReader( load_path )
%DATAREADER Allows to read the h5 data created by ImageIntegrationTool.
%   This function requires only the h5 file path. As an output it returns
%   vector of q values, time delays, ds_av matrix containing the averaged
%   differential curves for each of the time delay, ds_err matrix
%   containing the standard errors for each (q,t) point.
%   In addition, it provides Data structure, which contains total and
%   differential individual curves. Data is organized as follows:
%   Data.total.s_raw - unscaled total integrated curves
%   Data.total.s     - scaled total curves
%   Data.total.time  - time of image creation. It is a number of
%   seconds since 1st of January 0000
%   Data.total.scan  - scan number
%   Data.total.delay - time delay between laser and x-ray
%   Data.diff.ds     - individual differential curves 
%   Data.diff.time   - time when the LASER-ON curve was created
%   Data.diff.scan   - scan number
%   Data.diff.delay  - time delay of each differential curve
%   Data.diff.ds_av  - the same as ds_av
%   Data.diff.ds_err - the same as ds_err
%   Data.diff.t      - unique time delays
%   Data.diff.idx_incl - indexes of curves included into the calculation of
%   average curves.
%   Data.q           - q vector

try
Data.total.s_raw = h5read(load_path,'/total/s_raw');
Data.total.s     = h5read(load_path,'/total/s');
Data.total.time  = h5read(load_path,'/total/time');
Data.total.delay = h5read(load_path,'/total/delay');
    try
    Data.total.scan = h5read(load_path,'/total/scan');
    catch
    Data.diff.scan  = [];
    disp(['File ' load_path ' does not contain field /total/scan']);
    end
catch
Data.total.s_raw = [];
Data.total.s     = [];
Data.total.time  = [];
Data.total.delay = [];
Data.diff.scan   = [];
disp(['File ' load_path ' does not contain total scattering data']);
end

try
Data.diff.ds     = h5read(load_path,'/diff/ds');
Data.diff.time   = h5read(load_path,'/diff/time');
    try
    Data.diff.scan   = h5read(load_path,'/diff/scan');
    catch
    Data.diff.scan = [];
    disp(['File ' load_path ' does not contain field /diff/scan']);
    end
Data.diff.delay  = h5read(load_path,'/diff/delay');
Data.diff.ds_av  = h5read(load_path,'/diff/ds_av');
Data.diff.ds_err = h5read(load_path,'/diff/ds_err');
Data.diff.t      = h5read(load_path,'/diff/t');
Data.diff.idx_incl = h5read(load_path,'/diff/idx_incl');
catch
Data.diff.ds     = [];
Data.diff.time   = [];
Data.diff.scan = [];
Data.diff.delay  = [];
Data.diff.ds_av  = [];
Data.diff.ds_err = [];
Data.diff.t      = [];
Data.diff.idx_incl = [];
disp(['File ' load_path ' does not contain diff scattering data']);
end

try
Data.q           = h5read(load_path,'/q');
catch
Data.q = [];
disp('no q field was found');
end

ds_av = Data.diff.ds_av;
ds_err = Data.diff.ds_err;
q = Data.q;
t = Data.diff.t;

try
idx_incl = Data.diff.idx_incl;
Data.diff.idx_incl = cell(length(t),1);
for ii = 1:length(t)
    Data.diff.idx_incl{ii} = idx_incl(ii, idx_incl(ii,:)~=0);
end
catch
disp('no inclusion indexes available');
end
    
end


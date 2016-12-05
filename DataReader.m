function [  q, t, ds_av, ds_err, Data] =  DataReader( load_path )
%DATAREADER Summary of this function goes here
%   Detailed explanation goes here


Data.total.s_raw = h5read(load_path,'/total/s_raw');
Data.total.s     = h5read(load_path,'/total/s');
Data.total.time  = h5read(load_path,'/total/time');
try
Data.total.scan  = h5read(load_path,'/total/scan');
catch
Data.diff.scan = [];
disp(['File ' load_path ' does not contain field /total/scan']);
end
Data.total.delay = h5read(load_path,'/total/delay');

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

Data.q           = h5read(load_path,'/q');

ds_av = Data.diff.ds_av;
ds_err = Data.diff.ds_err;
q = Data.q;
t = Data.diff.t;

idx_incl = Data.diff.idx_incl;
Data.diff.idx_incl = cell(length(t),1);

for ii = 1:length(t)
    Data.diff.idx_incl{ii} = idx_incl(ii, idx_incl(ii,:)~=0);
end

end


function [ mask ] = readMask( old_mask_file )
    
idx_ext = strfind(old_mask_file,'.');
extension = old_mask_file(idx_ext+1:end);

if strcmp(extension,'h5')
    mask = h5read(old_mask_file,'/mask');
elseif strcmp(extension,'mccd')
    mask = imread(old_mask_file);
elseif strcmp(extension,'npy')
    mask = readNPY(old_mask_file);
end
mask = logical(mask);

end


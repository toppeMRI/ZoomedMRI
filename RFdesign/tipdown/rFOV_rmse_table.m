xy_fov = [24, 12, 10, 8];
z_under = [1, 2, 3, 6];

for ixy = 1:length(xy_fov)
    for iz = 1:length(z_under)
        fname = sprintf('rFOV_%d_x-y_%d_z',xy_fov(ixy),z_under(iz));
        load(fname); 
        RMSEs(ixy, iz) = rmse;  
    end
end
RMSEs

for ixy = 1:length(xy_fov)
    for iz = 1:length(z_under)
        fname = sprintf('rFOV_spgr_%d_x-y_%d_z',xy_fov(ixy),z_under(iz));
        load(fname); 
        RMSEs_spgr(ixy, iz) = rmse;  
    end
end
RMSEs_spgr

for ixy = 1:length(xy_fov)
    for iz = 1:length(z_under)
        fname = sprintf('rFOV_bssfp_%d_x-y_%d_z',xy_fov(ixy),z_under(iz));
        load(fname); 
        RMSEs_bssfp(ixy, iz) = rmse;  
    end
end
RMSEs_bssfp


for ixy = 1:length(xy_fov)
    for iz = 1:length(z_under)
        eval(['timeOneLeaf(ixy, iz) = size(', sprintf('k%dall',xy_fov(ixy)), ',1)*4e-3']) ; %ms 
        eval(['time(ixy, iz) = prod(size(', sprintf('k%dall',xy_fov(ixy)), '))*42/z_under(iz)*4e-3']) ; %ms 
    end
end
% timeOneLeaf
% time
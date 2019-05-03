function [ st_fn_shear_DIM ] = fn_ww__ext__dim__create_dim_shear( st_fn_shear_x, st_fn_shear_y )
%fn_ww__ext__dim__create_dim_shear: EXTERNAL Convert shear profile format
% 
%   [ st_fn_shear_DIM ] = fn_ww__ext__dim__create_dim_shear( st_fn_shear_x, st_fn_shear_y )
% 
% Convert our shear profile format to something we can use with DIM
% 
% TAGS: EXT
%
% See also
%   fn_ww__ext__dim__c()

st_fn_shear_DIM = struct;

st_fn_shear_DIM.fn_Ux = @(z) st_fn_shear_x.fn_phy_U(z);
st_fn_shear_DIM.fn_dUx = @(z) st_fn_shear_x.fn_phy_dU(z);
st_fn_shear_DIM.fn_ddUx = @(z) st_fn_shear_x.fn_phy_ddU(z);

st_fn_shear_DIM.fn_Uy = @(z) st_fn_shear_y.fn_phy_U(z);
st_fn_shear_DIM.fn_dUy = @(z) st_fn_shear_y.fn_phy_dU(z);
st_fn_shear_DIM.fn_ddUy = @(z) st_fn_shear_y.fn_phy_ddU(z);


end
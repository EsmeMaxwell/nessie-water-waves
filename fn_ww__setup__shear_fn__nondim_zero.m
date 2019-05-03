function [ st_fn_shear_zero, st_p ] = fn_ww__setup__shear_fn__nondim_zero( st_p, phy_h )
%fn_ww__setup__shear_fn__nondim_zero: Setup shear profile that is zero everywhere
%
%   [ st_fn_shear_zero, st_p ] = fn_ww__setup__shear_fn__nondim_zero( st_p, phy_h )
%
% Setup zero shear profile
%
% INPUT
%   st_p : Parameter set
%   phy_h : Dimensional depth
%
% OUTPUT
%   st_fn_shear_zero : Struct of shear profile functions
%   st_p : Update parameter set
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
% 
% See also
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__shear_fn__nondim_cospwr(),
%   fn_ww__setup__shear_fn__nondim_exp(),
%   fn_ww__setup__shear_fn__nondim_columbia_poly(),
%   fn_ww__setup__shear_fn__nondim_powerlaw()



phy_U0 = 0;

% Apply the nondimensionalisation to the parameter set
[ st_p ] = fn_ww__setup__param_ctl__re_cl( st_p, phy_U0, phy_h );

st_fn_shear_zero.fn_U = @(z) 0*z;
st_fn_shear_zero.fn_dU = @(z) 0*z;                            
st_fn_shear_zero.fn_ddU = @(z) 0*z;
                                
st_fn_shear_zero.fn_phy_U = @(z) 0*z;
st_fn_shear_zero.fn_phy_dU = @(z) 0*z;
st_fn_shear_zero.fn_phy_ddU = @(z)  0*z;
                               

% st_fn_shear_pwr.U_min = 0;
% st_fn_shear_pwr.U_max = 0;
% st_fn_shear_pwr.z_min = 0;
% st_fn_shear_pwr.z_max = 0;
% 
% st_fn_shear_cospwr.z_phy_min = 0;
% st_fn_shear_cospwr.z_phy_max = 0;
% st_fn_shear_cospwr.U_phy_min = 0;
% st_fn_shear_cospwr.U_phy_max = 0;

end
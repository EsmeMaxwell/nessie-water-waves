function [ st_fn_shear_exp, st_p ] = fn_ww__setup__shear_fn__nondim_exp( st_p, Fr2, h, phy_U0 )
%fn_ww__setup__shear_fn__nondim_exp: Setup shear profile for exponential profile
%
%   [ st_fn_shear_exp, st_p ] = fn_ww__setup__shear_fn__nondim_exp( st_p, Fr2, h, phy_U0 )
%
% Setup shear exponential shear profile. This is different from the other
% shear profiles in that we use d as a characteristic length scale.
%
% NOTE / WARNING : This has seriously annoying implications, must take
% phy_h suitably large but not too large as to mess up our calculations.
%
% INPUT
%   st_p : Parameter set
%   Fr2 : Desired shear Froude number
%   h : Depth
%   phy_U0 : Background current surface velocity
%   
% OUTPUT
%   st_fn_shear_exp : Struct of shear profile functions
%   st_p : Update parameter set
%
% TAGS: CORE, WWERRINSHEAR
% 
% See also
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__shear_fn__nondim_cospwr(),
%   fn_ww__setup__shear_fn__nondim_columbia_poly(),
%   fn_ww__setup__shear_fn__nondim_powerlaw(),
%   fn_ww__setup__shear_fn__nondim_zero()



assert( Fr2 > 0 && Fr2 <= 0.5, 'Check argument order, may have Fr2 too large.' );

% Apply the nondimensionalisation to the parameter set (this time we are
% also using phy_h as nondim h).
[ st_p ] = fn_ww__setup__param_ctl__re_cl__exp( st_p, Fr2, h, phy_U0 );

% Create actual shear function (nondimensional)
st_fn_shear_exp.fn_U = @(z) exp( z );
st_fn_shear_exp.fn_dU = @(z) exp( z );
st_fn_shear_exp.fn_ddU = @(z) exp( z );

% Do physical version
st_fn_shear_exp.fn_phy_U = @(z) phy_U0 * exp( z / st_p.phy_d );
st_fn_shear_exp.fn_phy_dU = @(z) phy_U0 * ( 1 / st_p.phy_d ) * exp( z / st_p.phy_d );
st_fn_shear_exp.fn_phy_ddU = @(z) phy_U0 * ( 1 / st_p.phy_d )^2 * exp( z / st_p.phy_d );


% % Calculate function U_max, U_min
% [ st_fn_shear_exp.U_min, st_fn_shear_exp.U_max, st_fn_shear_exp.z_min, st_fn_shear_exp.z_max ] = fn_ww__util__find_U_min_max( st_fn_shear_exp.fn_U, st_p );
% st_fn_shear_exp.z_phy_min = st_fn_shear_exp.z_min * st_p.phy_d;
% st_fn_shear_exp.z_phy_max = st_fn_shear_exp.z_max * st_p.phy_d;
% st_fn_shear_exp.U_phy_min = st_fn_shear_exp.fn_phy_U( st_fn_shear_exp.z_phy_min );
% st_fn_shear_exp.U_phy_max = st_fn_shear_exp.fn_phy_U( st_fn_shear_exp.z_phy_max );


end
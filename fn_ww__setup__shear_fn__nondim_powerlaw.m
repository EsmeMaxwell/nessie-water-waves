function [ st_fn_shear_pwr, st_p ] = fn_ww__setup__shear_fn__nondim_powerlaw( st_p, phy_U0, phy_h )
%fn_ww__setup__shear_fn__nondim_powerlaw: Setup shear profile for Colombia River
%
%   [ st_fn_shear_pwr, st_p ] = fn_ww__setup__shear_fn__nondim_powerlaw( st_p, phy_U0, phy_h )
%
% Setup 1/7th powerlaw shear profile
%
% INPUT
%   st_p : Parameter set
%   phy_U0 :  Background current surface velocity
%   phy_h : Dimensional depth
%
% OUTPUT
%   st_fn_shear_pwr : Struct of shear profile functions
%   st_p : Update parameter set
%
% TAGS: CORE, WWERRINSHEAR
% 
% See also
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__shear_fn__nondim_cospwr(),
%   fn_ww__setup__shear_fn__nondim_exp(),
%   fn_ww__setup__shear_fn__nondim_columbia_poly(),
%   fn_ww__setup__shear_fn__nondim_zero()




% Select a small bump factor because otherwise the bottom of the range has
% undefined derivatives
alpha = 1 - 1e-3;

% Apply the nondimensionalisation to the parameter set
[ st_p ] = fn_ww__setup__param_ctl__re_cl( st_p, phy_U0, phy_h );

% Create actual shear function (nondimensional). We write this natively
% nondimensional coordinates with U0=h=1 because the physical coordinates
% only define Fr2 in this instance.
st_fn_shear_pwr = struct;
st_fn_shear_pwr.fn_U = @(z) real( (alpha*z+1).^(1/7) );
st_fn_shear_pwr.fn_dU = @(z) real( (alpha/7) * (alpha*z+1).^(-6/7) );
st_fn_shear_pwr.fn_ddU = @(z) real( -(6*alpha^2/49) * (alpha*z+1).^(-13/7) );


st_fn_shear_pwr.fn_phy_U = @(z) real( phy_U0 * (alpha * ( z / phy_h ) + 1 ).^(1/7) );
st_fn_shear_pwr.fn_phy_dU = @(z) real( ( 1 / phy_h ) * phy_U0 * (alpha/7) * (alpha * ( z / phy_h ) + 1 ).^(-6/7) );
st_fn_shear_pwr.fn_phy_ddU = @(z) real( ( 1 / phy_h )^2 * phy_U0 * -(6*alpha^2/49) * (alpha * ( z / phy_h ) + 1 ).^(-13/7) );

% % Calculate function U_max, U_min
% [ st_fn_shear_pwr.U_min, st_fn_shear_pwr.U_max, st_fn_shear_pwr.z_min, st_fn_shear_pwr.z_max ] = fn_ww__util__find_U_min_max( st_fn_shear_pwr.fn_U, st_p );
% st_fn_shear_pwr.z_phy_min = st_fn_shear_pwr.z_min * phy_h;
% st_fn_shear_pwr.z_phy_max = st_fn_shear_pwr.z_max * phy_h;
% st_fn_shear_pwr.U_phy_min = st_fn_shear_pwr.fn_phy_U( st_fn_shear_pwr.z_phy_min );
% st_fn_shear_pwr.U_phy_max = st_fn_shear_pwr.fn_phy_U( st_fn_shear_pwr.z_phy_max );

% % Sanity check
% st_fn_shear_pwr.fn_U(0) - 1
% st_fn_shear_pwr.fn_dU(0) - 1/7
% 
% st_fn_shear_pwr.fn_dU(-1)
% st_fn_shear_pwr.fn_ddU(-1)
% error()

end
function [ st_fn_shear_poly, st_p ] = fn_ww__setup__shear_fn__nondim_columbia_poly( st_p )
%fn_ww__setup__shear_fn__nondim_columbia_poly: Setup shear profile for Colombia River
%
%   [ st_fn_shear_poly, st_p ] = fn_ww__setup__shear_fn__nondim_columbia_poly( st_p )
%
% Setup shear Colombia River shear profile
%
% INPUT
%   st_p : Parameter set
%
% OUTPUT
%   st_fn_shear_poly : Struct of shear profile functions
%   st_p : Update parameter set
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
% 
% See also
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__shear_fn__nondim_cospwr(),
%   fn_ww__setup__shear_fn__nondim_exp(),
%   fn_ww__setup__shear_fn__nondim_powerlaw(),
%   fn_ww__setup__shear_fn__nondim_zero()



phy_h = 20;
st_p.phy_h = phy_h;

% Raw polynomial (measurement coordinates)
v_cr_poly = -1 * [ -0.000003187151638 -0.000196584671643  -0.004694223054598  -0.054764717940654  -0.325720404172964   -0.978118228113837  -1.728173746026418 ];

% Readjust bottom value to be zero (this is assumed in s2.2 of the
% err-in-shear paper)
v_phys_poly = v_cr_poly;
v_phys_poly(end) = v_phys_poly(end) - polyval( v_phys_poly, - phy_h );
warning( 'Check why we zero this' );

% Maximum of polynomial (we must assume it's at the surface for now)
phy_U0 = polyval( v_phys_poly, 0 );

% Rescale
target_U0 = sqrt( 0.05 * st_p.phy_g * phy_h );
rescale_value = target_U0 / phy_U0;
v_phys_poly = v_phys_poly * rescale_value;

% Recalculate U0
phy_U0 = polyval( v_phys_poly, 0 );

% Create derivatives of physical poly
v_phys_d_poly = polyder( v_phys_poly );
v_phys_dd_poly = polyder( v_phys_d_poly );

% Must now convert to nondimensional coordinates (must rescale on vertical
% and use chainrule to do for z substitution). This is because we have
% $z = \tilde{z}h$ and so when converting form phy to nondim by
% substitution in the poly, we end up with the powers-of-h.
v_poly = ( 1 / phy_U0 ) * v_phys_poly .* [ phy_h^6 phy_h^5 phy_h^4 phy_h^3 phy_h^2 phy_h 1 ];
v_d_poly = polyder( v_poly );
v_dd_poly = polyder( v_d_poly );
% v_d_poly = ( 1 / phy_U0 ) * v_phys_d_poly .* [ phy_h^5 phy_h^4 phy_h^3 phy_h^2 phy_h 1 ];
% v_dd_poly = ( 1 / phy_U0 ) * v_phys_dd_poly .* [ phy_h^4 phy_h^3 phy_h^2 phy_h 1 ];

% Apply the nondimensionalisation to the parameter set
[ st_p ] = fn_ww__setup__param_ctl__re_cl( st_p, phy_U0, phy_h );

% Create actual shear function (nondimensional)
st_fn_shear_poly.fn_U = @(z) polyval( v_poly, z );
st_fn_shear_poly.fn_dU = @(z) polyval( v_d_poly, z );
st_fn_shear_poly.fn_ddU = @(z) polyval( v_dd_poly, z );

% Create shear profile (physical)
st_fn_shear_poly.fn_phy_U = @(z) polyval( v_phys_poly, z );
st_fn_shear_poly.fn_phy_dU = @(z) polyval( v_phys_d_poly, z );
st_fn_shear_poly.fn_phy_ddU = @(z) polyval( v_phys_dd_poly, z );
% st_fn_shear_poly.fn_phy_dU = @(z) ( 1 / st_p.phy_h ) * polyval( v_phys_d_poly, z );
% st_fn_shear_poly.fn_phy_ddU = @(z) ( 1 / st_p.phy_h )^2 * polyval( v_phys_dd_poly, z );


% % Calculate function U_max, U_min
% [ st_fn_shear_poly.U_min, st_fn_shear_poly.U_max, st_fn_shear_poly.z_min, st_fn_shear_poly.z_max ] = fn_ww__util__find_U_min_max( st_fn_shear_poly.fn_U, st_p );
% st_fn_shear_poly.z_phy_min = st_fn_shear_poly.z_min * phy_h;
% st_fn_shear_poly.z_phy_max = st_fn_shear_poly.z_max * phy_h;
% st_fn_shear_poly.U_phy_min = st_fn_shear_poly.fn_phy_U( st_fn_shear_poly.z_phy_min );
% st_fn_shear_poly.U_phy_max = st_fn_shear_poly.fn_phy_U( st_fn_shear_poly.z_phy_max );


end
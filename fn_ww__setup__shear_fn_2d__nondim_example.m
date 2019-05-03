function [ st_fn_shear_x, st_fn_shear_y, st_p ] = fn_ww__setup__shear_fn_2d__nondim_example( st_p )
%FN_WW__SETUP__SHEAR_FN_2D__NONDIM_EXAMPLE: Shear profile (2d example)
%
%   [ st_fn_shear_x, st_fn_shear_y, st_p ] = FN_WW__SETUP__SHEAR_FN_2D__NONDIM_EXAMPLE( st_p )
%

% Mental note: phy_U0 represents the supremum of the shear, in any
% direction; we use it as the characteristic velocity scale.  The Fr2
% number is calculated from that.  We can, hopefully, just scale using
% that.  To make the y direction weaker, use the beta scaling in our shear
% profile.


% Colombia river on x, cospwr on y


[ st_p ] = fn_ww__setup__param_std__re_cl(  );
phy_h = 20; Fr2 = 0.05;

[ st_fn_shear_cospwr, st_p ] = fn_ww__setup__shear_fn__nondim_cospwr( st_p, sqrt( Fr2 * st_p.phy_g * phy_h ), phy_h, 2, 4 * pi, 1, 0.5 );
[ st_fn_shear_poly, st_p ] = fn_ww__setup__shear_fn__nondim_columbia_poly( st_p );

st_fn_shear_x = st_fn_shear_cospwr;
st_fn_shear_y = st_fn_shear_poly;





% phy_h = 10;
% Fr2 = 0.005;
% phy_U0 = sqrt( Fr2 * st_p_ctl.phy_g * phy_h );
% 
% 
% alpha_x = 1.5;
% beta_x = 1.0;
% alpha_y = 0.6;
% beta_y = 0.5;
% 
% 
% 
% 
% % Apply the nondimensionalisation to the parameter set, assume we're using
% % sqrcos on the x axis and that it is stronger
% [ st_p_ctl ] = fn_ww__setup__param_ctl__re_cl( phy_U0, phy_h );
% [ st_p ] = fn_ww__setup__merge_parameters( st_p, st_p_ctl );
% 
% 

% % Sqrcos on the x axis alpha = 1.5*pi
% st_fn_shear_x.fn_U = @(z) beta_x * (1/2) * ( cos( 2* alpha_x * z ) + 1 );
% st_fn_shear_x.fn_dU = @(z) -alpha_x * beta_x * sin( 2 * alpha_x * z );
% st_fn_shear_x.fn_ddU = @(z) -2 * alpha_x.^2 * beta_x * cos( 2 * alpha_x * z );
% 
% st_fn_shear_x.fn_phy_U = @(z) phy_U0 * beta_x * (1/2) * ( cos( 2 * alpha_x * z / phy_h ) + 1 );
% st_fn_shear_x.fn_phy_dU = @(z) ( 1 / phy_h ) * phy_U0 * -alpha_x * beta_x * sin( 2 * alpha_x * z / phy_h );
% st_fn_shear_x.fn_phy_ddU = @(z) ( 1 / phy_h )^2 * phy_U0 * -2 * alpha_x.^2 * beta_x * cos( 2 * alpha_x * z / phy_h );
% 
% 
% % Sqrcos on y but with lower alpha and lower shear
% st_fn_shear_y.fn_U = @(z) beta_y * (1/2) * ( cos( 2* alpha_y * z ) + 1 );
% st_fn_shear_y.fn_dU = @(z) -alpha_y * beta_y * sin( 2 * alpha_y * z );
% st_fn_shear_y.fn_ddU = @(z) -2 * alpha_y.^2 * beta_y * cos( 2 * alpha_y * z );
% 
% st_fn_shear_y.fn_phy_U = @(z) phy_U0 * beta_y * (1/2) * ( cos( 2 * alpha_y * z / phy_h ) + 1 );
% st_fn_shear_y.fn_phy_dU = @(z) ( 1 / phy_h ) * phy_U0 * -alpha_y * beta_y * sin( 2 * alpha_y * z / phy_h );
% st_fn_shear_y.fn_phy_ddU = @(z) ( 1 / phy_h )^2 * phy_U0 * -2 * alpha_y.^2 * beta_y * cos( 2 * alpha_y * z / phy_h );




end
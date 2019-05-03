function [ st_fn_shear_cospwr, st_p ] = fn_ww__setup__shear_fn__nondim_cospwr( st_p, phy_U0, phy_h, alpha, beta, gamma, delta )
%fn_ww__setup__shear_fn__nondim_cospwr: Setup shear profile U_T (for cos of z to power)
%
%   [ st_fn_shear_cospwr, st_p ] = fn_ww__setup__shear_fn__nondim_cospwr( st_p, phy_U0, phy_h, alpha, beta, gamma, delta )
%
% Setup shear profile U_T.
%
% INPUT
%   st_p : Parameter set
%   phy_U0 : Dimensional background current surface velocity
%   phy_h : Dimsional depth
%   alpha, beta, gamma, delta : Parameters for the shear profile
%
% OUTPUT
%   st_fn_shear_cospwr : Struct of shear profile functions
%   st_p : Update parameter set
%
% TAGS: CORE, SISCPFLIB
% 
% EXAMPLE
%   [ st_p ] = fn_ww__setup__param_std__re_cl(  );
%   phy_h = 20; Fr2 = 0.05;
%   [ st_fn_shear, st_p ] = fn_ww__setup__shear_fn__nondim_cospwr( st_p, sqrt( Fr2 * st_p.phy_g * phy_h ), phy_h, 2, 4 * pi, 1, 0.5 );
% 
% See also
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__shear_fn__nondim_columbia_poly(),
%   fn_ww__setup__shear_fn__nondim_exp(),
%   fn_ww__setup__shear_fn__nondim_powerlaw(),
%   fn_ww__setup__shear_fn__nondim_zero()



% Apply the nondimensionalisation to the parameter set
[ st_p ] = fn_ww__setup__param_ctl__re_cl( st_p, phy_U0, phy_h );

st_fn_shear_cospwr.fn_U = @(z) 0.5 * gamma * ( 1 + delta * z ) .* cos( beta * (-z).^alpha ) + 0.5;

st_fn_shear_cospwr.fn_dU = @(z) 0.5 * gamma * ( delta * cos( beta * (-z).^alpha ) ...
                                + alpha * beta * (-z) .^ ( alpha - 1 ) .* ( 1 + delta * z ) .* sin( beta * (-z).^alpha ) );
                            
st_fn_shear_cospwr.fn_ddU = @(z) 0.5 * gamma * ( 2 * alpha * beta * delta * (-z) .^ ( alpha - 1 ) .* sin( beta * (-z).^alpha ) ...
                                    + ( 1 + delta * z ) .* ( -alpha^2 * beta^2 * (-z) .^ ( 2 * ( alpha - 1 ) ) .* cos( beta * (-z).^alpha ) ...
                                    -alpha * ( alpha - 1 ) * beta * (-z) .^ ( alpha - 2 ) .* sin( beta * (-z).^alpha ) ) );


st_fn_shear_cospwr.fn_phy_U = @(z) 0.5 * phy_U0 * gamma * ( 1 + delta * ( z / phy_h ) ) .* cos( beta * ( -z / phy_h ).^alpha ) + 0.5 * phy_U0;

st_fn_shear_cospwr.fn_phy_dU = @(z) 0.5 * ( 1 / phy_h ) * phy_U0 * gamma * ( ...
                                delta * cos( beta * ( -z / phy_h ).^alpha ) ...
                                + alpha * beta * ( -z / phy_h ) .^ ( alpha - 1 ) .* ( 1 + delta * z / phy_h ) .* sin( beta * ( -z / phy_h ).^alpha ) ...
                                            );
                            
st_fn_shear_cospwr.fn_phy_ddU = @(z) 0.5 * ( 1 / phy_h )^2 * phy_U0 * gamma * ( 2 * alpha * beta * delta * ( -z / phy_h ) .^ ( alpha - 1 ) .* sin( beta * ( -z / phy_h ).^alpha ) ...
                                    + ( 1 + delta * z / phy_h ) .* ( -alpha^2 * beta^2 * ( -z / phy_h ) .^ ( 2 * ( alpha - 1 ) ) .* cos( beta * ( -z / phy_h ).^alpha ) ...
                                    -alpha * ( alpha - 1 ) * beta * ( -z / phy_h ) .^ ( alpha - 2 ) .* sin( beta * ( -z / phy_h ).^alpha ) ) );


                                
% % Calculate function U_max, U_min
% [ st_fn_shear_cospwr.U_min, st_fn_shear_cospwr.U_max, st_fn_shear_cospwr.z_min, st_fn_shear_cospwr.z_max ] = fn_ww__util__find_U_min_max( st_fn_shear_cospwr.fn_U, st_p );
% st_fn_shear_cospwr.z_phy_min = st_fn_shear_cospwr.z_min * phy_h;
% st_fn_shear_cospwr.z_phy_max = st_fn_shear_cospwr.z_max * phy_h;
% st_fn_shear_cospwr.U_phy_min = st_fn_shear_cospwr.fn_phy_U( st_fn_shear_cospwr.z_phy_min );
% st_fn_shear_cospwr.U_phy_max = st_fn_shear_cospwr.fn_phy_U( st_fn_shear_cospwr.z_phy_max );

                                
end
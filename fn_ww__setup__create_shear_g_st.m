function [ st_g_shear ] = fn_ww__setup__create_shear_g_st( st_fn_shear_x, st_fn_shear_y, v_theta, st_p )
%fn_ww__setup__create_shear_r: Setup return a reduced shear profile struct
%
%   [ st_r_shear ] = fn_ww__setup__create_shear_r( st_fn_shear )
%
% Almost trivial, just stores the shear fns in another struct along with
% Umin, Umax results.
%
% INPUT
%   st_fn_shear_x, st_fn_shear_y : Structs of shear profile functions
%
% OUTPUT
%   st_r_shear : Struct containing shear functions and Umin,Umax inf
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__shear_fn__nondim_cospwr(),
%   fn_ww__setup__shear_fn__nondim_exp(),
%   fn_ww__setup__shear_fn__nondim_columbia_poly(),
%   fn_ww__setup__shear_fn__nondim_powerlaw(),
%   fn_ww__setup__shear_fn__nondim_zero()


% Copy shear functions
st_g_shear = struct;
st_g_shear.st_fn_shear_x = st_fn_shear_x;
st_g_shear.st_fn_shear_y = st_fn_shear_y;

% Store v_theta as we're appropriate only to this v_theta now
st_g_shear.v_theta = v_theta;

% Calculate function U_max, U_min
Ntheta = numel( v_theta );
ca_r_shear = cell( 1, Ntheta );

for lp_theta=1:Ntheta
    
    theta = v_theta( lp_theta );
    
    st_fn_shear = struct;
    st_fn_shear.fn_U = @(z) cos(theta) * st_fn_shear_x.fn_U(z) + sin(theta) * st_fn_shear_y.fn_U(z);
    st_fn_shear.fn_dU = @(z) cos(theta) * st_fn_shear_x.fn_dU(z) + sin(theta) * st_fn_shear_y.fn_dU(z);
    st_fn_shear.fn_ddU = @(z) cos(theta) * st_fn_shear_x.fn_ddU(z) + sin(theta) * st_fn_shear_y.fn_ddU(z);
    st_fn_shear.fn_phy_U = @(z) cos(theta) * st_fn_shear_x.fn_phy_U(z) + sin(theta) * st_fn_shear_y.fn_phy_U(z);
    st_fn_shear.fn_phy_dU = @(z) cos(theta) * st_fn_shear_x.fn_phy_dU(z) + sin(theta) * st_fn_shear_y.fn_phy_dU(z);
    st_fn_shear.fn_phy_ddU = @(z) cos(theta) * st_fn_shear_x.fn_phy_ddU(z) + sin(theta) * st_fn_shear_y.fn_phy_ddU(z);
    
    [ st_crit ] = fn_ww__util__find_U_min_max( st_fn_shear, st_p );
        
    ca_r_shear{lp_theta}.st_fn_shear = st_fn_shear;
    ca_r_shear{lp_theta}.st_crit = st_crit;
    
end


st_g_shear.ca_r_shear = ca_r_shear;

end
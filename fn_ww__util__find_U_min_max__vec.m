function [ st_crit ] = fn_ww__util__find_U_min_max__vec( st_vec, st_v_shear, st_p )
%fn_ww__util__find_U_min_max: Util find Umin and Umax of vector shear profile
% 
%   [ U_min, U_max ] = fn_ww__util__find_U_min_max( v_U )
% 
% 
% 

st_crit = struct;

% Find min, max
[ st_crit.U_min, st_crit.U_min_idx ] = min( st_v_shear.v_U );
[ st_crit.U_max, st_crit.U_max_idx ] = max( st_v_shear.v_U );

% Now identify the z coordinates for this
st_crit.z_min = st_vec.v_zm( st_crit.U_min_idx );
st_crit.z_max = st_vec.v_zm( st_crit.U_max_idx );

% Do for physical coords if required
if ( isfield( st_v_shear, 'v_phy_U' ) )
    st_crit.U_phy_min = st_v_shear.v_phy_U( st_crit.U_min_idx );
    st_crit.U_phy_max = st_v_shear.v_phy_U( st_crit.U_max_idx );
    st_crit.z_phy_min = st_vec.v_zp( st_crit.U_min_idx );
    st_crit.z_phy_max = st_vec.v_zp( st_crit.U_max_idx );
end


end
function [ st_v_shear ] = fn_ww__setup__shear_fn_to_vec( st_vec, st_fn_shear, st_p )
%fn_ww__setup__shear_fn_to_vec: Setup discretise continuous profile to vector
%
%   [ st_v_shear ] = fn_ww__setup__shear_fn_to_vec( st_Dn, st_fn_shear, st_p )
%
% Take a vector of points and supplied shear profile and returns vector
% form of shear profile, sampled at those points. Specifically, it will
% contain vectors v_U, v_dU, and v_ddU.
%
% INPUT
%   st_vec : Struct of vector parts of st_Dn, can be st_Dn if wanted.
%   st_fn_shear : Struct of shear profile functions
%
% OUTPUT
%   st_v_shear : Struct of vector versions of the shear profile
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__setup__diffmtrx__WR_chebdif(),
%   fn_ww__setup__diffmtrx__WR_poldif(),
%   fn_ww__setup__diffmtrx_mp__WR_chebdif(),
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__shear_fn__nondim_cospwr(),
%   fn_ww__setup__shear_fn__nondim_exp(),
%   fn_ww__setup__shear_fn__nondim_columbia_poly(),
%   fn_ww__setup__shear_fn__nondim_powerlaw(),
%   fn_ww__setup__shear_fn__nondim_zero()



% Nondimensional first
st_v_shear = struct;
st_v_shear.v_U = st_fn_shear.fn_U( st_vec.v_zm );
st_v_shear.v_dU = st_fn_shear.fn_dU( st_vec.v_zm );
st_v_shear.v_ddU = st_fn_shear.fn_ddU( st_vec.v_zm );


% Physical coordinates, we only do this if physical vector exists
if ( isfield( st_vec, 'v_zp' ) & isfield( st_fn_shear, 'fn_phy_U' ) )

    st_v_shear.v_phy_U = st_fn_shear.fn_phy_U( st_vec.v_zp );
    st_v_shear.v_phy_dU = st_fn_shear.fn_phy_dU( st_vec.v_zp );
    st_v_shear.v_phy_ddU = st_fn_shear.fn_phy_ddU( st_vec.v_zp );

end



end
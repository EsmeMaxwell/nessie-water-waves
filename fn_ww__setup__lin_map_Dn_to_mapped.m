function [ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p )
%fn_ww__setup__lin_map_Dn_to_mapped: Setup map diff mtrcs+vecs to nondimesional+dimensional coordinates
%
%   [ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, a, b )
%
% Creates the simple linearly mapped differentiation matrices and z vector.
%
%
% INPUT
%   st_Dn : Standard diff matrices and vector in [-1,1] domain.
%   st_p : Parameters
% 
% OUTPUT
%   st_Dn : Same struct but with physical fields filled in. 
% 
% 
% See also
%   fn_ww__setup__lin_conv_z0_to_zm(),
%   fn_ww__setup__diffmtrx__WR_chebdif(),
%   fn_ww__setup__diffmtrx_mp__WR_chebdif(),
%   fn_ww__setup__diffmtrx__WR_poldif


% Nondim coordinates
st_Dn.v_zm = fn_ww__setup__lin_conv_z0_to_zm( st_Dn.v_z0, st_p.a, st_p.b );
st_Dn.a_D2m = ( 4 / ( st_p.b - st_p.a )^2 ) * st_Dn.a_D2;
st_Dn.a_Dm = ( 2 / ( st_p.b - st_p.a ) ) * st_Dn.a_D;

% Physical coordinates
st_Dn.v_zp = fn_ww__setup__lin_conv_z0_to_zm( st_Dn.v_z0, st_p.phy_a, st_p.phy_b );
st_Dn.a_D2p = ( 4 / ( st_p.phy_b - st_p.phy_a )^2 ) * st_Dn.a_D2;
st_Dn.a_Dp = ( 2 / ( st_p.phy_b - st_p.phy_a ) ) * st_Dn.a_D;



end
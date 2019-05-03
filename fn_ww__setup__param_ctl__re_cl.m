function [ st_p ] = fn_ww__setup__param_ctl__re_cl( st_p, phy_U0, phy_h )
%fn_ww__setup__param_ctl__re_cl: Setup finish parameter setup with known shear fn with nondim against fixed depth
%
%   [ st_p ] = fn_ww__setup__param_ctl__re_cl( phy_U0, phy_h )
% 
% Finishes setup of parameter set and nondimensionalisation, etc, after
% shear profile is known. This alters how the physical properties are
% setup, how mappings will work, and Fr2.
%
% For similar function that uses d as characteristic length for exp
% profile, see fn_ww__setup__param_ctl__re_cl__exp().
%
% This should generally only be called within shear setup functions like
% fn_ww__setup__shear_fn__nondim_cospwr()
%
% NOTE : This assumes the characteristic length scale is h!
%
% INPUT
%   phy_U0 :  Background current surface velocity
%   phy_h : Dimensional depth
%
% OUTPUT
%   st_p_ctl : Parameters merged into the main parameter set
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__setup__shear_fn__nondim_cospwr(),
%   fn_ww__setup__param_ctl__re_cl__exp(),
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__shear_fn_to_vec()


% These are used in a merge to overwrite
st_p_ctl.h = 1;
st_p_ctl.phy_h = phy_h;
st_p_ctl.a = -st_p_ctl.h;
st_p_ctl.b = 0;
st_p_ctl.phy_a = -phy_h;
st_p_ctl.phy_b = 0;

st_p_ctl.phy_U0 = phy_U0;
st_p_ctl.Fr2 = st_p_ctl.phy_U0.^2 / ( st_p.phy_g  * st_p_ctl.phy_h );

% Save characteristic lengthscale
st_p_ctl.charlen = st_p_ctl.phy_h;

% Need to adjust safety threshold up
st_p_ctl.fp_safety_thrsh = 200000;

% % Ensure PF initialistation is done in MP
% st_p_ctl.bp_pf_use_eig_mp = true;

[ st_p ] = fn_ww__setup__merge_parameters( st_p, st_p_ctl );


end
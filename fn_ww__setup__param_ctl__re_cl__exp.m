function [ st_p ] = fn_ww__setup__param_ctl__re_cl__exp( st_p, Fr2, h, phy_U0 )
%fn_ww__setup__param_ctl__re_cl__exp: Return standard parametersa RE CL with exp profile
%
%   [ st_p ] = fn_ww__setup__param_ctl__re_cl__exp( Fr2, h, phy_U0 )
%
% Finishes setup of parameter set and nondimensionalisation, etc, after
% shear profile is known. This alters how the physical properties are
% setup, how mappings will work, and Fr2.
%
% NOTES : Assumes characeristic length scale is d (related to exponential
% profile) and that h is in some sense "large".
%
% For similar function that uses h as characteristic length for exp
% profile, see fn_ww__setup__param_ctl__re_cl().
%
% This should generally only be called within shear setup functions like
% fn_ww__sim_shear_prtb__full_exp_profile()
%
% INPUT
%   Fr2 : Shear Froude number to use
%   h : Depth
%   phy_U0 :  Background current surface velocity
%
% OUTPUT
%   st_p_ctl : Parameters to be merged into the main parameter set
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__sim_shear_prtb__full_exp_profile(),
%   fn_ww__setup__param_ctl__re_cl__exp(),
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__shear_fn_to_vec()

% We use these early on... (this should really be in one place but have
% coded self into corener, fix later)
st_p_ctl.phy_g = 9.80665;

% Nondim setup first, as we based everything else off this
st_p_ctl.h = h;
st_p_ctl.a = -st_p_ctl.h;
st_p_ctl.b = 0;
st_p_ctl.Fr2 = Fr2;

% Physical domain setup
st_p_ctl.phy_U0 = phy_U0;
st_p_ctl.phy_d = st_p_ctl.phy_U0^2 / ( st_p_ctl.phy_g * Fr2 );
st_p_ctl.phy_h = h * st_p_ctl.phy_d;
st_p_ctl.phy_a = -st_p_ctl.phy_h;
st_p_ctl.phy_b = 0;

% Save characteristic lengthscale
st_p_ctl.charlen = st_p_ctl.phy_d;

% Need to adjust safety threshold up
st_p_ctl.fp_safety_thrsh = 200000;

% % Ensure PF initialistation is done in MP
% st_p_ctl.bp_pf_use_eig_mp = true;

[ st_p ] = fn_ww__setup__merge_parameters( st_p, st_p_ctl );


end
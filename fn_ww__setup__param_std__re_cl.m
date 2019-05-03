function [ st_p ] = fn_ww__setup__param_std__re_cl(  )
%fn_ww__setup__param_std__re_cl: Return standard parameters for RE CL method
%
%   [ st_p ] = fn_ww__setup__param_std__re_cl(  )
% 
% Setup code to create initial parameter set, should usually be done first.
% Subsequent setup should usually be done indirectly, calling one of the
% shear setup functions, e.g. fn_ww__setup__shear_fn__nondim_cospwr(),
% rather than calling fn_ww__setup__param_ctl__re_cl() or
% fn_ww__setup__param_ctl__re_cl__exp() directly.
%
% OUTPUT
%   st_p : Struct of default parameters
%
% See also,
%   fn_ww__setup__param_ctl__re_cl(),
%   fn_ww__setup__param_ctl__re_cl__exp(),
%   fn_ww__setup__shear_fn__nondim_cospwr(),
%   fn_ww__setup__shear_fn__nondim_exp(),
%   fn_ww__setup__shear_fn_to_vec()


% We use these early on... 
st_p.phy_g = 9.80665;


% For now, we're going to force the user to generate the basic parameters
% using the shear profile specific parameter creation... avoids potential
% accidents.


% Some tricks to handle a few numerical issues
st_p.fp_safety_thrsh = 200;
st_p.fp_imag_thrsh = 1e-14;
st_p.fp_current_shift = 0;
%st_p.bp_eigenvector_filter = true;
%st_p.bp_delete_zero_eig = false;
st_p.fp_zero_tol = 1e-14;
st_p.bp_fld = false;
st_p.bp_mp = false;
st_p.bp_disp_update = true;
st_p.ip_mp_digits = 34;
st_p.ip_mp_threads = 0;
st_p.bp_phy_calc = false;
st_p.bp_ex_btm_bc = true;   % Whether to remove bottom Dirichlet boundary condition (marginally better results)
st_p.ip_std_evec_norm = 100;    % Use surface velocity as standard eigenvector normalisation
st_p.bp_pf_force_unit_evec = true;
st_p.bp_pf_debug = false;

st_p.bp_Uminmax_vec = false;   % Whether to use vector calculations for vector Umin,Umax calculations... it's faster but less stable.
st_p.bp_Uminmax_precalc = true; % Whether to precalc, true should force st_p.bp_Uminmax_vec to false!!
if ( st_p.bp_Uminmax_precalc )
    st_p.bp_Uminmax_vec = false;
end

st_p.fp_adtpv_depth_tol = eps;  % Tolerance at which we consider the eigenfunction to be effectively zero
st_p.fp_adtpv_depth_hmin = 0.3;     % Proportionate target minimum depth
st_p.fp_adtpv_depth_hmax = 0.8;     % Proportionate target maximum depth
st_p.ip_adtpv_min_numel_overlap = 5;    % Minimum number of elements required in overlap region
st_p.fp_adtpv_depth_tau = 0.98;

st_p.bp_apx_force_vec_shear = false;    % Whether to force use of vector shear profile in approximation routines.

end
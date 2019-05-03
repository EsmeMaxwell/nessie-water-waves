function [ st_fit_olscb ] = fn_ww__calc_fit__ols( st_fit_cb_interp, st_fit_olscb_mtrx, v_zs_candidate_U, h )
%fn_ww__calc_fit__ols: Calc OLS fit in Chebyshev basis
%
%   [ st_fit_olscb ] = fn_ww__calc_fit__ols( st_fit_cb_interp, st_fit_olscb_mtrx, v_zs_candidate_U, h )
%
% Calc OLS using Chebyshev polynomial basis. Note that the st_cb_interp
% basis is for the original domain (it was implicitly mapped to the
% Chebyshev domain in fn_ww__calc_fit__prep_cb_basis(), we just undo that
% here).
%
% TAGS: WWERRINSHEAR
%
% See also
%   fn_ww__calc_fit__prep_cb_basis(),
%   fn_ww__calc_fit__prep_ols_matrices(),
%   fn_ww__calc_fit__prep_lin_surf_extrapolate()


% Note: need the cb basis to include z=0 but not use it in the original
% fit.


% Calculate the mean and std dev of the sample and temporarily rescale and
% translate
mean_candidate_U = mean( v_zs_candidate_U );
std_candidate_U = std( v_zs_candidate_U );
v_zs_candidate_rescale = ( v_zs_candidate_U - mean_candidate_U ) / std_candidate_U;

v_beta = st_fit_olscb_mtrx.a_V * st_fit_olscb_mtrx.a_Spinv * (st_fit_olscb_mtrx.a_U') * v_zs_candidate_rescale;

v_zs_fit_U_rescale = st_fit_cb_interp.a_T * v_beta;
v_zs_fit_dU_rescale = ( 2 / h ) * st_fit_cb_interp.a_dT * v_beta;
v_zs_fit_ddU_rescale = ( 4 / h^2 ) * st_fit_cb_interp.a_ddT * v_beta;


% Rescale and translate back to original before returning
v_zs_fit_U = ( v_zs_fit_U_rescale * std_candidate_U ) + mean_candidate_U;
v_zs_fit_dU = ( v_zs_fit_dU_rescale * std_candidate_U );
v_zs_fit_ddU = ( v_zs_fit_ddU_rescale * std_candidate_U );

st_fit_olscb = struct;
st_fit_olscb.v_zs_fit_U = v_zs_fit_U;   % Review March 2019: this is the shear profile over the required Chebyshev pts
st_fit_olscb.v_zs_fit_dU = v_zs_fit_dU;
st_fit_olscb.v_zs_fit_ddU = v_zs_fit_ddU;
st_fit_olscb.v_beta = v_beta;



end
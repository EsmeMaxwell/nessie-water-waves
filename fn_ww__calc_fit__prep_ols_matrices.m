function [ st_fit_olscb_mtrx ] = fn_ww__calc_fit__prep_ols_matrices( st_fit_cb )
%fn_ww__calc_fit__prep_ols_matrices: Calc OLS Chebyshev precalc matrices
%
%   [ st_fit_olscb_mtrx ] = fn_ww__calc_fit__prep_ols_matrices( st_fit_cb )
%
% From CB basis, calculates necessary SVD decomposition for use in ordinary
% least squares
%
% TAGS: WWERRINSHEAR
%
% See also
%   fn_ww__calc_fit__prep_cb_basis(),
%   fn_ww__calc_fit__ols(),
%   fn_ww__calc_fit__prep_lin_surf_extrapolate()



st_fit_olscb_mtrx = struct;

% For OLS using SVD
[ a_U, a_S, a_V ] = svd( st_fit_cb.a_T );
a_Spinv = ( a_S )';
a_Spinv( a_Spinv ~= 0 ) = 1 ./ a_Spinv( a_Spinv ~= 0 );  % Eugh! Fix me.

st_fit_olscb_mtrx.a_U = a_U;
st_fit_olscb_mtrx.a_S = a_S;
st_fit_olscb_mtrx.a_V = a_V;
st_fit_olscb_mtrx.a_Spinv = a_Spinv;


% v_beta_2 = a_V * a_Spinv * (a_U') * v_zm_sample_rand_nd;
% v_zm_sample_ols_cb_2_nl = a_A * v_beta_2;


end
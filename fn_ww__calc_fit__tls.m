function [ st_fit_tlscb ] = fn_ww__calc_fit__tls( st_fit_cb_interp, st_fit_tlscb_mtrx, v_zs_candidate_U, h )
%fn_ww__calc_fit__tls: Calc TLS fit in Chebyshev basis
%
%   [ st_fit_tlscb ] = fn_ww__calc_fit__tls( st_fit_cb_interp, st_fit_tlscb_mtrx, v_zs_candidate_U, h )
%
% Calc TLS using Chebyshev polynomial basis. Note that the st_cb_interp
% basis is for the original domain (it was implicitly mapped to the
% Chebyshev domain in fn_ww__calc_fit__prep_cb_basis(), we just undo that
% here).
%
% TAGS: WWERRINSHEAR
% 
% See also
%   fn_ww__calc_fit__prep_cb_basis(),
%   fn_ww__calc_fit__prep_tls_matrices(),
%   fn_ww__calc_fit__prep_lin_surf_extrapolate()


% Calculate the mean and std dev of the sample and temporarily rescale and
% translate
mean_candidate_U = mean( v_zs_candidate_U );
std_candidate_U = std( v_zs_candidate_U );
v_zs_candidate_rescale = ( v_zs_candidate_U - mean_candidate_U ) / std_candidate_U;

% Do the calculation
v_zs_candidate_rescale_Tk = [ v_zs_candidate_rescale; zeros( st_fit_tlscb_mtrx.rows_T, 1 ) ];
v_beta = st_fit_tlscb_mtrx.a_V * st_fit_tlscb_mtrx.a_Spinv * (st_fit_tlscb_mtrx.a_U') * v_zs_candidate_rescale_Tk;

v_zs_fit_U_rescale = st_fit_cb_interp.a_T * v_beta;
v_zs_fit_dU_rescale = ( 2 / h ) * st_fit_cb_interp.a_dT * v_beta;
v_zs_fit_ddU_rescale = ( 4 / h^2 ) * st_fit_cb_interp.a_ddT * v_beta;


% Rescale and translate back to original before returning
v_zs_fit_U = ( v_zs_fit_U_rescale * std_candidate_U ) + mean_candidate_U;
v_zs_fit_dU = ( v_zs_fit_dU_rescale * std_candidate_U );
v_zs_fit_ddU = ( v_zs_fit_ddU_rescale * std_candidate_U );

st_fit_tlscb = struct;
st_fit_tlscb.v_zs_fit_U = v_zs_fit_U;
st_fit_tlscb.v_zs_fit_dU = v_zs_fit_dU;
st_fit_tlscb.v_zs_fit_ddU = v_zs_fit_ddU;
st_fit_tlscb.v_beta = v_beta;



end
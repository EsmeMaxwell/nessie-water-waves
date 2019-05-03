function [ v_zs_candidate_U_extr ] = fn_ww__calc_fit__prep_lin_surf_extrapolate( v_zs, v_zs_candidate_U, Npts )
%fn_ww__calc_fit__prep_lin_surf_extrapolate: Calc shear profile extrapolation from near surface to surface
% 
%   [ v_zs_candidate_U_extr ] = fn_ww__calc_fit__prep_lin_surf_extrapolate( v_zs, v_zs_candidate_U, Npts )
%
% Assuming v_zs does not reach the surface, we use linear extrapolation of
% the Npts topmost points to obtain a guess of the surface point. The
% reconstructed v_zs and candidate is returned, which now includes the
% estimate of the surface current. This means that polynomial fitting can
% then be done and we don't risk massive overshoot or cause a discontinuous
% derivative.
%
% TAGS: WWERRINSHEAR
%
%

% Basic sanity checks
assert( isnumeric( v_zs ) && 1 == size( v_zs, 2 ), 'Column vector expected' );
assert( numel( v_zs ) > Npts, 'Not enough points to extrapolate with' );
assert( v_zs(1) < 0, 'Vecter already has a surface entry' );

% Use MATLAB polyfit to handle this
st_poly = polyfit( v_zs(1:Npts), v_zs_candidate_U(1:Npts), 1 );
surf_val = polyval( st_poly, [ 0 ] );

%v_zs_extr = [ 0; v_zs ];
v_zs_candidate_U_extr = [ surf_val; v_zs_candidate_U ];


end
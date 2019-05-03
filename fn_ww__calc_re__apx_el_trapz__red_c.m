function [ v_c_EL_trapz ] = fn_ww__calc_re__apx_el_trapz__red_c( v_k, Nz, st_fn_shear, st_p )
%fn_ww__calc_re__apx_el_trapz__red_c: Calc reduced DR approximation Ellingsen-Li using trapz quadrature
% 
%   [ v_c_EL_trapz ] = fn_ww__calc_re__apx_el_trapz__red_c( v_k, Nz, st_fn_shear, st_p )
% 
% Calculates an approximation to the dispersion relation for Rayeligh
% equation with free-surface using Ellingsen-Li approximation scheme in the
% reduced problem setting.  The integral is evaluated using trapz. You may
% prefer to use fn_ww__calc_re__apx_el_cc__red_c().
%
% TAGS: SISCPFLIB
% 
% INPUT
%
%   v_k : k vector
%   Nz : Number of z points used
%   st_fn_shear : Shear profile
%   st_p : Parameter set
%
% See also
%   fn_ww__calc_re__apx_el_cc__red_c()


Nk = numel( v_k );

v_z = linspace( st_p.a, st_p.b, Nz ).';

v_c0 = sqrt( ( 1./v_k ) .* tanh( st_p.h * v_k ) );

a_z = repmat( v_z, 1, Nk );
a_dUz = repmat( st_fn_shear.fn_dU( v_z ), 1, Nk );

a_delta_integrand = a_dUz .* ( exp( 2 * v_k .* a_z ) - exp( -4 * st_p.h * v_k - 2 * v_k .* a_z ) ) ./ ( 1 - exp( -4 * st_p.h * v_k ) );
v_delta = trapz( v_z, a_delta_integrand ) ./ v_c0;


v_c1 = v_c0 .* ( sqrt( 1 + v_delta.^2 ) - v_delta );

v_c_EL_trapz = st_fn_shear.fn_U( 0 ) + v_c1;

%v_c_EL2_trapz =  - ( ( v_c0 .* v_c1 ) ./ ( 2 * sqrt( 1 + v_delta.^2 ) ) );


end
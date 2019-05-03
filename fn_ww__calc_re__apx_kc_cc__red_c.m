function [ v_c_KC_cc ] = fn_ww__calc_re__apx_kc_cc__red_c( v_k, v_z_cc, v_zw_cc, st_fn_shear, st_p )
%fn_ww__calc_re__apx_kc_cc__red_c: Calc reduced DR approximation Kirby-Chen using Curtis-Clenshaw quadrature
% 
%   [ v_c_KC_cc ] = fn_ww__calc_re__apx_kc_cc__red_c( v_k, v_z_cc, v_zw_cc, st_fn_shear, st_p )
% 
% Calculates an approximation to the dispersion relation for Rayeligh
% equation with free-surface using Kirby-Chen approximation scheme in the
% reduced problem setting.  The integral is evaluated using high order
% Curtis-Clenshaw quadrature.
%
% TAGS: SISCPFLIB
% 
% INPUT
%
%   v_k : k vector
%   v_z_cc : CC absiccas (calculated from fn_ww__ext__cc_weights())
%   v_zw_cc : CC weights (calculated from fn_ww__ext__cc_weights())
%   st_fn_shear : Shear profile
%   st_p : Parameter set
%
% See also
%   fn_ww__ext__cc_weights()


Nk = numel( v_k );
Nz = numel( v_z_cc );

if ( st_p.bp_phy_calc )
    v_c0 = sqrt( ( st_p.phy_g ./ v_k ) .* tanh( st_p.phy_h * v_k ) );

    a_z_cc = repmat( v_z_cc, 1, Nk );
    a_zw_cc = repmat( v_zw_cc, 1, Nk );
    a_Uz_cc = repmat( st_fn_shear.fn_phy_U( v_z_cc ), 1, Nk );

    a_c1_integrand = a_Uz_cc .* ( exp( 2 * v_k .* a_z_cc ) + exp( -4 * st_p.phy_h * v_k - 2 * v_k .* a_z_cc ) ) ./ ( 1 - exp( -4 * st_p.phy_h * v_k ) );
    v_c1_ccint = sum( a_zw_cc .* a_c1_integrand );
    v_c1 = 2 * v_k .* v_c1_ccint;

    v_c_KC_cc = v_c0 + v_c1;
else
    v_c0 = sqrt( ( 1 ./ ( v_k * st_p.Fr2 ) ) .* tanh( st_p.h * v_k ) );

    a_z_cc = repmat( v_z_cc, 1, Nk );
    a_zw_cc = repmat( v_zw_cc, 1, Nk );
    a_Uz_cc = repmat( st_fn_shear.fn_U( v_z_cc ), 1, Nk );

    a_c1_integrand = a_Uz_cc .* ( exp( 2 * v_k .* a_z_cc ) + exp( -4 * st_p.h * v_k - 2 * v_k .* a_z_cc ) ) ./ ( 1 - exp( -4 * st_p.h * v_k ) );
    v_c1_ccint = sum( a_zw_cc .* a_c1_integrand );
    v_c1 = 2 * v_k .* v_c1_ccint;

    v_c_KC_cc = v_c0 + v_c1;    
end

end
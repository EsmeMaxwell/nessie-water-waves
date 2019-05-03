function [ v_c_EL_cc ] = fn_ww__calc_re__apx_el_cc__red_c( v_k, v_z_cc, v_zw_cc, st_r_shear, st_p )
%fn_ww__calc_re__apx_el_cc__red_c: Calc reduced DR approximation Ellingsen-Li using Curtis-Clenshaw quadrature
% 
%   [ v_c_EL_cc ] = fn_ww__calc_re__apx_el_cc__red_c( v_k, v_z_cc, v_zw_cc, st_fn_shear, st_p )
% 
% Calculates an approximation to the dispersion relation for Rayleigh
% equation with free-surface using Ellingsen-Li approximation scheme in the
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
    
    % Deal with function or vector shear profile stuff first
    if ( isfield( st_r_shear, 'st_fn_shear' ) & ~st_p.bp_apx_force_vec_shear )
        v_shear_phy_dU_cc = st_r_shear.st_fn_shear.fn_phy_dU( v_z_cc );
    else
        assert( numel( v_z_cc ) == numel( st_r_shear.st_v_Shear.v_phy_dU ), 'Incorrect vector lengths!' );
        v_shear_phy_dU_cc = st_r_shear.st_v_shear.v_phy_dU;
    end
         
    % Do the calculations...
    v_c0 = sqrt( ( st_p.phy_g ./ v_k ) .* tanh( st_p.phy_h * v_k ) );

    a_z_cc = repmat( v_z_cc, 1, Nk );
    a_zw_cc = repmat( v_zw_cc, 1, Nk );
    a_dUz_cc = repmat( v_shear_phy_dU_cc, 1, Nk );

    a_delta_integrand = a_dUz_cc .* ( exp( 2 * v_k .* a_z_cc ) - exp( -4 * st_p.phy_h * v_k - 2 * v_k .* a_z_cc ) ) ./ ( 1 - exp( -4 * st_p.phy_h * v_k ) );
    v_delta = sum( a_zw_cc .* a_delta_integrand ) ./ v_c0;

    v_c1 = v_c0 .* ( sqrt( 1 + v_delta.^2 ) - v_delta );

    v_c_EL_cc = st_fn_shear.fn_phy_U( 0 ) + v_c1;
    
else
    
    % Deal with function or vector shear profile stuff first
    if ( isfield( st_r_shear, 'st_fn_shear' ) & ~st_p.bp_apx_force_vec_shear )
        v_shear_dU_cc = st_r_shear.st_fn_shear.fn_dU( v_z_cc );
        shear_U0 = st_r_shear.st_fn_shear.fn_U( 0 );
    else
        assert( numel( v_z_cc ) == numel( st_r_shear.st_v_shear.v_dU ), 'Incorrect vector lengths!' );
        v_shear_dU_cc = st_r_shear.st_v_shear.v_dU;
        shear_U0 = st_r_shear.st_v_shear.v_U(1);
    end    
    
    % Do the calculations...
    v_c0 = sqrt( ( 1 ./ ( v_k * st_p.Fr2 ) ) .* tanh( st_p.h * v_k ) );

    a_z_cc = repmat( v_z_cc, 1, Nk );
    a_zw_cc = repmat( v_zw_cc, 1, Nk );
    a_dUz_cc = repmat( v_shear_dU_cc, 1, Nk );

    a_delta_integrand = a_dUz_cc .* ( exp( 2 * v_k .* a_z_cc ) - exp( -4 * st_p.h * v_k - 2 * v_k .* a_z_cc ) ) ./ ( 1 - exp( -4 * st_p.h * v_k ) );
    v_delta = sum( a_zw_cc .* a_delta_integrand ) ./ v_c0;

    v_c1 = v_c0 .* ( sqrt( 1 + v_delta.^2 ) - v_delta );

    v_c_EL_cc = shear_U0 + v_c1;

end



end
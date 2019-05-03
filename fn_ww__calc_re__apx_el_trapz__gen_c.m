function [ v_c_EL_trapz ] = fn_ww__calc_re__apx_el_trapz__gen_c( v_k, v_theta, Nz, st_fn_shear_x, st_fn_shear_y, st_p )
%fn_ww__calc_re__apx_el_trapz__gen_c: DEV
%

Nk = numel( v_k );

v_z = linspace( st_p.a, st_p.b, Nz ).';

[ v_kx, v_ky ] = pol2cart( v_theta, v_k );


v_c0 = sqrt( ( 1./v_k ) .* tanh( st_p.h * v_k ) );

a_z = repmat( v_z, 1, Nk );
a_dUx_z = repmat( st_fn_shear_x.fn_dU( v_z ), 1, Nk );
a_dUy_z = repmat( st_fn_shear_y.fn_dU( v_z ), 1, Nk );

a_delta_integrand = ( v_kx .* a_dUx_z + v_ky .* a_dUy_z ).* ( exp( 2 * v_k .* a_z ) - exp( -4 * st_p.h * v_k - 2 * v_k .* a_z ) ) ./ ( 1 - exp( -4 * st_p.h * v_k ) );
v_delta = trapz( v_z, a_delta_integrand ) ./ v_c0;


v_c1 = v_c0 .* ( sqrt( 1 + v_delta.^2 ) - v_delta );

v_c_EL_trapz = ( v_kx * st_fn_shear_x.fn_U( 0 ) + v_ky * st_fn_shear_y.fn_U( 0 ) ) + v_c1;

%v_c_EL2_trapz =  - ( ( v_c0 .* v_c1 ) ./ ( 2 * sqrt( 1 + v_delta.^2 ) ) );


end
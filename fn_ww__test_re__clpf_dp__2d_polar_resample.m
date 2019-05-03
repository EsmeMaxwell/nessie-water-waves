function [ st_gen_resample ] = fn_ww__test_re__clpf_dp__2d_polar_resample( st_p )
%fn_ww__test_re__clpf_dp__2d_polar_resample: DEV


k_min = 0.005;
k_max = 80;
k0 = 2.0;
tol = 1e-12;

Nz = 48;
Nthetaresample = 20;
Nkresample = 30;

Nthetaq = 100;
Nkq = 100;

[ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( Nz, 1 );
[ st_Dn_mp ] = fn_ww__setup__diffmtrx_mp__WR_chebdif( Nz );

[ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p.a, st_p.b );
[ st_Dn_mp ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_mp, st_p.a, st_p.b );

[ st_fn_shear_x ] = fn_ww__setup__shear_fn( 10, { 'csp_b', 0.05 } );
[ st_fn_shear_y ] = fn_ww__setup__shear_fn( 10, { 'csp_b', 0.025 } );

[ st_v_shear_x ] = fn_ww__setup__shear_fn_to_vec( st_Dn.v_zm, st_fn_shear_x );
[ st_v_shear_y ] = fn_ww__setup__shear_fn_to_vec( st_Dn.v_zm, st_fn_shear_y );


[ st_gen_resample ] = fn_ww__calc_re__clpf_dp__gen_2d_polar_resample_c( st_Dn, k_min, k_max, Nkresample, Nthetaresample, st_v_shear_x, st_v_shear_y, k0, tol, st_p );


% Create polar query points
v_theta_q = fn_ww__util__create_theta_vec( Nthetaq );
v_k_q = fn_ww__util__create_k_vec( k_min, k_max, Nkq, 3, 1e-4 );
[ v_theta_q_rs, v_k_q_rs, v_kx_q_rs, v_ky_q_rs ] = fn_ww__util__create_grid_ordered_pair_polar( v_k_q, v_theta_q, true );

[ v_c_PF_interp ] = fn_ww__calc_re__clpf_interp__2d_resample_c( st_gen_resample, v_theta_q_rs, v_k_q_rs );


figure(1);
scatter3( v_kx_q_rs, v_ky_q_rs, v_c_PF_interp, 'bo' );


end
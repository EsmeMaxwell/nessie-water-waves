function [ st_res ] = fn_ww__test_re__cl__2d_basic(  )
%fn_ww__test_re__cl__2d_basic: DEV

k_min = 0.025;
k_max = 250;
Nk = 60;
Ntheta = 10;
Nz = 64;


[ v_k ] = fn_ww__util__create_k_vec( k_min, k_max, Nk, 1, 0 );
[ v_theta ] = fn_ww__util__create_theta_vec( Ntheta );



[ st_p ] = fn_ww__setup__param_std__re_cl(  );
[ st_fn_shear_x, st_fn_shear_y, st_p ] = fn_ww__setup__shear_fn_2d__nondim_example( st_p )

[ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( Nz, 1 );
[ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p );

[ st_v_shear_x ] = fn_ww__setup__shear_fn_to_vec( st_Dn, st_fn_shear_x, st_p );
[ st_v_shear_y ] = fn_ww__setup__shear_fn_to_vec( st_Dn, st_fn_shear_y, st_p );

[ st_fn_shear_DIM ] = fn_ww__ext__dim__create_dim_shear( st_fn_shear_x, st_fn_shear_y );

a_c_rad = zeros( Ntheta, Nk );
a_c_ang = zeros( Ntheta, Nk );

% Do radial first
for lp_theta=1:Ntheta
    
    fprintf( '\ntheta %d of %d\n', lp_theta, Ntheta );
    
    [ a_c_rad( lp_theta, : ) ] = fn_ww__calc_re__cl__gen_sl_pol_rad_c( st_Dn, v_k, v_theta(lp_theta), st_v_shear_x, st_v_shear_y, st_p );    
end



% Do angular
for lp_k=1:Nk
    
    fprintf( '\nk %d of %d\n', lp_k, Nk);
    
    [ v_c_ang ] = fn_ww__calc_re__cl__gen_sl_pol_ang_c( st_Dn, v_k(lp_k), v_theta, st_v_shear_x, st_v_shear_y, st_p );
    a_c_ang( :, lp_k ) = v_c_ang.';
end



% can probably use fn_ww__util__create_grid_ordered_pair_polar() ?

[ a_k, a_theta ] = meshgrid( v_k, v_theta );
v_theta_rs = reshape( a_theta, 1, [] );
v_k_rs = reshape( a_k, 1, [] );
[ v_kx_rs, v_ky_rs ] = pol2cart( v_theta_rs, v_k_rs );
v_c_rad_rs = reshape( a_c_rad, 1, [] ); 
v_c_ang_rs = reshape( a_c_ang, 1, [] ); 



figure(1);
scatter3( v_kx_rs, v_ky_rs, v_c_rad_rs, 'bo' );
hold on;
scatter3( v_kx_rs, v_ky_rs, v_c_ang_rs, 'mx' );
hold off;

st_res = struct;


end
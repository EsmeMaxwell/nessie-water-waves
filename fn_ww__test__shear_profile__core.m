function [ st_res ] = fn_ww__test__shear_profile__core( st_r_shear, v_k, k_min, k_max, st_p, st_p_mp )
%fn_ww__test__shear_profile__core: Test some std uses of shear profile
%
%   [ st_res ] = fn_ww__test__shear_profile__core( st_fn_shear, v_k, k_min, k_max, st_p, st_p_mp )
%
% Called from other functions testing specific shear profiles, this
% function does generic tests on an arbitrary shear profile with provided
% parameters.
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__test__shear_profile__do_powerlaw(),
%   fn_ww__test__shear_profile__do_exp()


pf_tol = 1e-11;
k0 = 3.0;

% Normal DIM
DIMl_tol = 1e-10;
DIMl_z = 24;
DIMl_itr = 10;

% HA DIM
DIMh_tol = 1e-12;
DIMh_z = 10000;
DIMh_itr = 20000;



[ st_Dn_48 ] = fn_ww__setup__diffmtrx__WR_poldif( 48, 1 );
[ st_Dn_48 ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_48, st_p );
[ st_Dn_128 ] = fn_ww__setup__diffmtrx__WR_poldif( 128, 1 );
[ st_Dn_128 ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_128, st_p );
[ st_Dn_mp ] = fn_ww__setup__diffmtrx_mp__WR_chebdif( 160 );
[ st_Dn_mp ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_mp, st_p_mp );

v_k_phy = v_k / st_p.phy_h;




%---- Nondim calcs

fprintf( 'Nondim...\n' );

% ADPTV versions
[ st_adptv ] = fn_ww__calc_re__adptv_interval_refine( v_k, st_p );
st_proc = struct;

% ADPTV collocation
st_proc.ip_method = 1;
[ v_c_p_adtpv_48 ] = fn_ww__calc_re__adptv__red_c( st_Dn_48, st_r_shear, st_adptv, st_proc, st_p );
[ v_c_p_adtpv_128 ] = fn_ww__calc_re__adptv__red_c( st_Dn_128, st_r_shear, st_adptv, st_proc, st_p );
[ v_c_p_adtpv_mp ] = fn_ww__calc_re__adptv__red_c( st_Dn_mp, st_r_shear, st_adptv, st_proc, st_p_mp );

% ADPTV PF
st_proc.ip_method = 2;
st_proc.pf_tol = pf_tol;
[ v_c_PF_adtpv_48 ] = fn_ww__calc_re__adptv__red_c( st_Dn_48, st_r_shear, st_adptv, st_proc, st_p );
[ v_c_PF_adtpv_128 ] = fn_ww__calc_re__adptv__red_c( st_Dn_128, st_r_shear, st_adptv, st_proc, st_p );
[ v_c_PF_adtpv_mp ] = fn_ww__calc_re__adptv__red_c( st_Dn_mp, st_r_shear, st_adptv, st_proc, st_p_mp );


% Approx first
fprintf( '... approx ... \n');
[ v_z_cc, v_zw_cc ]=fn_ww__ext__cc_weights( 256, st_p.a, st_p.b );
[ v_c_EL_cc ] = fn_ww__calc_re__apx_el_cc__red_c( v_k, v_z_cc, v_zw_cc, st_r_shear.st_fn_shear, st_p );
[ v_c_KC_cc ] = fn_ww__calc_re__apx_kc_cc__red_c( v_k, v_z_cc, v_zw_cc, st_r_shear.st_fn_shear, st_p );

% Collocation
fprintf( '... CL ... \n');
[ v_c_p_48 ] = fn_ww__calc_re__cl__red_c( st_Dn_48, v_k, st_r_shear, st_p );
[ v_c_p_128 ] = fn_ww__calc_re__cl__red_c( st_Dn_128, v_k, st_r_shear, st_p );
[ v_c_p_mp ] = fn_ww__calc_re__cl__red_c( st_Dn_mp, v_k, st_r_shear, st_p_mp );

% PF
fprintf( '... PF ... \n');
[ st_dp_48 ] = fn_ww__calc_re__clpf_dp__red_c( st_Dn_48, k_min, k_max, st_r_shear, k0, pf_tol, st_p );
[ v_c_PF_48 ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_48, v_k, false );
[ st_dp_128 ] = fn_ww__calc_re__clpf_dp__red_c( st_Dn_128, k_min, k_max, st_r_shear, k0, pf_tol, st_p );
[ v_c_PF_128 ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_128, v_k, false );
[ st_dp_mp ] = fn_ww__calc_re__clpf_dp__red_c( st_Dn_mp, k_min, k_max, st_r_shear, k0, pf_tol, st_p_mp );
[ v_c_PF_mp ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_mp, v_k, false );



%---- Physical
st_p.bp_phy_calc = true;
fprintf( '\nPhysical...\n' )

% Approx first
fprintf( '... approx ... \n');
[ v_z_cc, v_zw_cc ]=fn_ww__ext__cc_weights( 256, st_p.phy_a, st_p.phy_b );
[ v_c_EL_cc_phy ] = fn_ww__calc_re__apx_el_cc__red_c( v_k_phy, v_z_cc, v_zw_cc, st_r_shear.st_fn_shear, st_p );
[ v_c_KC_cc_phy ] = fn_ww__calc_re__apx_kc_cc__red_c( v_k_phy, v_z_cc, v_zw_cc, st_r_shear.st_fn_shear, st_p );

% Basic collocation
fprintf( '... CL ... \n');
[ v_c_p_48_phy ] = fn_ww__calc_re__cl__red_c( st_Dn_48, v_k_phy, st_r_shear, st_p );

% DIM
fprintf( '... DIM ... \n');
st_fn_shear_DIM = struct;
st_fn_shear_DIM.fn_Ux = @(z) st_r_shear.st_fn_shear.fn_phy_U(z);
st_fn_shear_DIM.fn_dUx = @(z) st_r_shear.st_fn_shear.fn_phy_dU(z);
st_fn_shear_DIM.fn_ddUx = @(z) st_r_shear.st_fn_shear.fn_phy_ddU(z);
st_fn_shear_DIM.fn_Uy = @(z) 0*z;
st_fn_shear_DIM.fn_dUy = @(z) 0*z;
st_fn_shear_DIM.fn_ddUy = @(z) 0*z;

[ tc_KCl, tc_ELl, tc_dimMl, itr_cnt_shortl, itr_cnt_longl ] = fn_ww__ext__dim__c( v_k_phy, 0, st_fn_shear_DIM, DIMh_tol, DIMh_z, DIMh_itr, st_p.phy_h, st_p.phy_g, false );
[ tc_KCh, tc_ELh, tc_dimMh, itr_cnt_shorth, itr_cnt_longh ] = fn_ww__ext__dim__c( v_k_phy, 0, st_fn_shear_DIM, DIMh_tol, DIMh_z, DIMh_itr, st_p.phy_h, st_p.phy_g, false );



st_res = struct;

st_res.v_k = v_k;

st_res.v_c_EL_cc = v_c_EL_cc;
st_res.v_c_KC_cc = v_c_KC_cc;

st_res.v_c_p_48 = v_c_p_48;
st_res.v_c_p_128 = v_c_p_128;
st_res.v_c_p_mp = v_c_p_mp;

st_res.v_c_PF_48 = v_c_PF_48;
st_res.v_c_PF_128 = v_c_PF_128;
st_res.v_c_PF_mp = v_c_PF_mp;

st_res.v_c_p_adtpv_48 = v_c_p_adtpv_48;
st_res.v_c_p_adtpv_128 = v_c_p_adtpv_128;
st_res.v_c_p_adtpv_mp = v_c_p_adtpv_mp;

st_res.v_c_PF_adtpv_48 = v_c_PF_adtpv_48;
st_res.v_c_PF_adtpv_128 = v_c_PF_adtpv_128;
st_res.v_c_PF_adtpv_mp = v_c_PF_adtpv_mp;



st_res.v_c_EL_cc_phy = v_c_EL_cc_phy / st_p.phy_U0;
st_res.v_c_KC_cc_phy = v_c_KC_cc_phy / st_p.phy_U0;

st_res.v_c_p_48_phy = v_c_p_48_phy / st_p.phy_U0;

st_res.v_c_DIMl = tc_dimMl / st_p.phy_U0;
st_res.v_c_KCDIMl = tc_KCl / st_p.phy_U0;
st_res.v_c_ELDIMl = tc_ELl / st_p.phy_U0;
st_res.itr_cnt_shortl = itr_cnt_shortl;
st_res.itr_cnt_longl = itr_cnt_longl;

st_res.v_c_DIMh = tc_dimMh / st_p.phy_U0;
st_res.v_c_KCDIMh = tc_KCh / st_p.phy_U0;
st_res.v_c_ELDIMh = tc_ELh / st_p.phy_U0;
st_res.itr_cnt_shorth = itr_cnt_shorth;
st_res.itr_cnt_longh = itr_cnt_longh;

st_res.st_p = st_p;
st_res.st_p_mp = st_p_mp;


end
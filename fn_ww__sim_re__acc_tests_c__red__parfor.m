function [ v_c_CL, v_c_PF, v_c_PFmp, v_c_ADPTV_CL, v_c_ADPTV_PF, v_c_ADPTV_PFmp, v_c_DIM, st_err_CL, st_err_PF ] = fn_ww__sim_re__acc_tests_c__red__parfor( Nz, st_r_shear, st_adptv, k_min, k_max, v_k, v_k_phy, st_p )
%fn_ww__sim_re__acc_tests_c__red__parfor: Simulation (parfor fn) CL+PF eig accuracy, backwards error, cond number tests
% 
%   [ v_c_CL, v_c_PF, v_c_PFmp, v_c_DIM, st_err_CL, st_err_PF ] = fn_ww__sim_re__acc_tests_c__red__parfor( Nz, st_fn_shear, k_min, k_max, v_k, v_k_phy, st_p )
% 
% The parfor function called from fn_ww__sim_re__acc_tests_c__red().
% 
% TAGS: SISCPFLIB
%
% See also
%   fn_ww__sim_re__acc_tests_c__red()



PF_epsilon = 1e-5;
PF_tol = 5e-15; % Should probably note in MS that running PF with rediculously tight tolerance
k0 = 2.0;
mp_digits = 34;


% Setup diff matrices
[ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( Nz, 1 );
[ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p );


[ st_Dn_mp ] = fn_ww__setup__diffmtrx_mp__WR_chebdif( Nz, mp_digits );
[ st_Dn_mp ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_mp, st_p );
[ st_p_PFmp ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_pf_use_eig_mp', true, 'st_Dn_mp', st_Dn_mp, 'ip_mp_digits', mp_digits ) );



% Setup required parameters
[ st_p_CL ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_disp_update', false ) );    
[ st_p_PF ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_disp_update', false, 'bp_pf_debug', true ) );

% Do CL calc
[ v_c_CL, ~, v_residual, v_berr, v_cond ] = fn_ww__calc_re__cl__red_c( st_Dn, v_k, st_r_shear, st_p_CL );
st_err_CL = struct;
st_err_CL.v_residual = v_residual;
st_err_CL.v_berr = v_berr;
st_err_CL.v_cond = v_cond;

% Do PF calc
[ st_dp ] = fn_ww__calc_re__clpf_dp__red_c( st_Dn, k_min-PF_epsilon, k_max+PF_epsilon, st_r_shear, k0, PF_tol, st_p_PF );
[ v_c_PF, ~ ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp, v_k, 0 );
st_err_PF = struct;
st_err_PF.v_berr = st_dp.v_berr;
st_err_PF.v_cond = st_dp.v_cond;

% Do PFmp calc
[ st_dp_mp ] = fn_ww__calc_re__clpf_dp__red_c( st_Dn_mp, k_min-PF_epsilon, k_max+PF_epsilon, st_r_shear_mp, k0, PF_tol, st_p_PFmp );
[ v_c_PFmp, ~ ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_mp, v_k, 0 );


% Do CL+ADPTV calc

% Do PF+ADPTV calc


% Do PFmp+ADPTV calc



% Do DIM calc
[ st_fn_shear_zero, st_p_DIM ] = fn_ww__setup__shear_fn__nondim_zero( st_p, st_p.phy_h );
[ st_fn_shear_DIM ] = fn_ww__ext__dim__create_dim_shear( st_r_shear.st_fn_shear, st_fn_shear_zero );
[ ~, ~, v_c_DIM, ~, ~ ] = fn_ww__ext__dim__c( v_k_phy, 0, st_fn_shear_DIM, 1e-12, Nz, 10, st_p.phy_h, st_p.phy_g, false );
v_c_DIM = v_c_DIM / st_p.phy_U0;




end
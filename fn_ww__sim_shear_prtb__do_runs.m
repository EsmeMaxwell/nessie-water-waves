function [ st_err ] = fn_ww__sim_shear_prtb__do_runs( st_fn_shear_exact, v_zs, v_zs_err_width, v_k, v_percentile, Nitr, Nz, st_p, s_filename, s_backup_filename )
%fn_ww__sim_shear_prtb__do_runs: Sim stochatic runs for err-in-shear problem
%
%   fn_ww__sim_shear_prtb__do_runs( st_fn_shear_exact, v_zs, v_zs_err_width, v_k, v_percentile, Nitr, Nz, st_p, s_filename )
% 
% Performs the multiple runs requires for the error-in-shear-profile
% problem. Calls fn_ww__sim_shear_prtb__do_single() repeatedly.
%
% Should always be called by appropriate management function, e.g.
% fn_ww__sim_shear_prtb__full_columbia_profile(),
% fn_ww__sim_shear_prtb__full_exp_profile(), or
% fn_ww__sim_shear_prtb__full_powerlaw_profile()
% 
% TAGS: WWERRINSHEAR
% 
% See also
%   fn_ww__sim_shear_prtb__do_single()
%   fn_ww__sim_shear_prtb__full_columbia_profile(),
%   fn_ww__sim_shear_prtb__full_exp_profile(),
%   fn_ww__sim_shear_prtb__full_powerlaw_profile()


if ( nargout > 0 )
    st_err = struct;
    st_err.b_ok = true;
end

% We can, almost, use this to test how strong near-surface curvature may
% affect things.  Really only pertinant for the Bezier fit.
Nz_cc = 4 * Nz;

Nk = numel( v_k );
k_min = min( v_k ) - 1e-4;
k_max = max( v_k ) + 1e-4;
pf_tol = 1e-9;
k0 = 2.0;

cb_poly_deg = 7;  % need to check what degree of poly this actually creates
ml_poly_deg = 7; 

b_do_extr = ( 0 ~= v_zs(1) );

% Setup shear struct
[ st_r_shear_exact ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear_exact, st_p );

% Setup differentiation matrix, etc.
[ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( Nz, 1 );
[ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p );

% Setup abscissa and weights for approximation scheme
[ v_z_cc, v_zw_cc ] = fn_ww__ext__cc_weights( Nz_cc, st_p.a, st_p.b );

% Setup adptv regions
[ st_adptv ] = fn_ww__calc_re__adptv_interval_refine( v_k, st_p );

% Check whether we need to initialise an MP diff matrix
if ( isfield( st_p, 'bp_pf_use_eig_mp' ) && st_p.bp_pf_use_eig_mp )
    [ st_p.st_Dn_mp ] = fn_ww__setup__diffmtrx_mp__WR_chebdif( Nz );
    [ st_p.st_Dn_mp ] = fn_ww__setup__lin_map_Dn_to_mapped( st_p.st_Dn_mp, st_p );
end

% Calc exact vector
v_zs_exact_U = st_fn_shear_exact.fn_U( v_zs );

% % Setup reference (why, do we still use this?)
% [ st_v_shear_exact ] = fn_ww__setup__shear_fn_to_vec( st_Dn, st_fn_shear_exact, st_p );

% figure(1); plot( st_Dn.v_zm, st_v_shear_exact.v_U, 'b' );
% figure(2); plot( st_Dn.v_zm, st_v_shear_exact.v_dU, 'g' );
% figure(3); plot( st_Dn.v_zm, st_v_shear_exact.v_ddU, 'm' );
% st_v_shear_exact.v_U
% st_v_shear_exact.v_dU
% st_v_shear_exact.v_ddU
% error()
% figure(1); plot( st_v_shear_exact.v_U, st_Dn.v_zm, 'b' );
% error()

% Calculate reference data
% [ st_dp_ref ] = fn_ww__calc_re__clpf_dp__red_c( st_Dn, k_min, k_max, st_v_shear_exact, k0, pf_tol, st_p ); 
% [ v_c_ref ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_ref, v_k, false );
st_proc = struct;
st_proc.ip_method = 2;
st_proc.pf_tol = pf_tol;
[ v_c_ref ] = fn_ww__calc_re__adptv__red_c( st_Dn, st_r_shear_exact, st_adptv, st_proc, st_p );
[ v_apx_ref ] = fn_ww__calc_re__apx_el_cc__red_c( v_k, v_z_cc, v_zw_cc, st_r_shear_exact, st_p );

% figure(1);
% semilogx( v_k, v_c_ref, 'b', v_k, v_apx_EL_ref, 'm--' );
% figure(2);
% loglog( v_k, abs(v_c_ref - v_apx_EL_ref), 'b' );
% error()

% Prepare fit resources, i.e. bases, fit matrices, etc.
[ st_fit_cb ] = fn_ww__calc_fit__prep_cb_basis( cb_poly_deg, v_zs );
if ( b_do_extr )
    v_zs_extr = [ 0; v_zs ];
    [ st_fit_cb_interp_PF ] = fn_ww__calc_fit__prep_cb_basis( cb_poly_deg, st_Dn.v_zm, -v_zs(1) );
    [ st_fit_cb_interp_APX ] = fn_ww__calc_fit__prep_cb_basis( cb_poly_deg, v_z_cc, -v_zs(1) );
    
    [ st_fit_cb_extr ] = fn_ww__calc_fit__prep_cb_basis( cb_poly_deg, v_zs_extr );
    [ st_fit_olscb_mtrx_extr ] = fn_ww__calc_fit__prep_ols_matrices( st_fit_cb_extr );
    [ st_fit_tlscb_mtrx_extr ] = fn_ww__calc_fit__prep_tls_matrices( st_fit_cb_extr, 2, 1.0 );
    [ st_fit_cb_interp_PF_extr ] = fn_ww__calc_fit__prep_cb_basis( cb_poly_deg, st_Dn.v_zm );
    [ st_fit_cb_interp_APX_extr ] = fn_ww__calc_fit__prep_cb_basis( cb_poly_deg, v_z_cc );
else
    [ st_fit_cb_interp_PF ] = fn_ww__calc_fit__prep_cb_basis( cb_poly_deg, st_Dn.v_zm );
    [ st_fit_cb_interp_APX ] = fn_ww__calc_fit__prep_cb_basis( cb_poly_deg, v_z_cc );
end
[ st_fit_olscb_mtrx ] = fn_ww__calc_fit__prep_ols_matrices( st_fit_cb );
[ st_fit_tlscb_mtrx ] = fn_ww__calc_fit__prep_tls_matrices( st_fit_cb, 2, 1 );
[ fn__exp_fit ] = @( v_param, v_x ) v_param(3) + v_param(1) * exp( v_param(2) * v_x );
[ fn__exp_fit_d ] = @( v_param, v_x ) v_param(1) * v_param(2) * exp( v_param(2) * v_x );
[ fn__exp_fit_dd ] = @( v_param, v_x ) v_param(1) * v_param(2)^2 * exp( v_param(2) * v_x );
% [ fn__exp_fit ] = @( v_param, v_x ) v_param(1) * exp( v_param(2) * v_x );
% [ fn__exp_fit_d ] = @( v_param, v_x ) v_param(1) * v_param(2) * exp( v_param(2) * v_x );
% [ fn__exp_fit_dd ] = @( v_param, v_x ) v_param(1) * v_param(2)^2 * exp( v_param(2) * v_x );


st_res = struct;
st_res.v_z_cc = v_z_cc;
st_res.v_zw_cc = v_zw_cc;
st_res.st_adptv = st_adptv;
st_res.cb_poly_deg = cb_poly_deg;
st_res.ml_poly_deg = ml_poly_deg;
st_res.st_fn_shear_exact = st_fn_shear_exact;
st_res.st_fit_cb = st_fit_cb;
st_res.st_fit_olscb_mtrx = st_fit_olscb_mtrx;
st_res.st_fit_tlscb_mtrx = st_fit_tlscb_mtrx;
st_res.st_fit_cb_interp_PF = st_fit_cb_interp_PF;
st_res.st_fit_cb_interp_APX = st_fit_cb_interp_APX;
st_res.k_min = k_min;
st_res.k_max = k_max;
st_res.pf_tol = pf_tol;
st_res.k0 = k0;
st_res.fn__exp_fit = fn__exp_fit;
st_res.fn__exp_fit_d = fn__exp_fit_d;
st_res.fn__exp_fit_dd = fn__exp_fit_dd;
st_res.v_c_ref = v_c_ref;
st_res.v_apx_ref = v_apx_ref;
if ( b_do_extr )
    st_res.st_fit_cb_extr = st_fit_cb_extr;
    st_res.st_fit_olscb_mtrx_extr = st_fit_olscb_mtrx_extr;
    st_res.st_fit_tlscb_mtrx_extr = st_fit_tlscb_mtrx_extr;
    st_res.v_zs_extr = v_zs_extr;
    st_res.st_fit_cb_interp_PF_extr = st_fit_cb_interp_PF_extr;
    st_res.st_fit_cb_interp_APX_extr = st_fit_cb_interp_APX_extr;
end


a_c_err_lsqexp_PF = zeros( Nitr, Nk, 'mp' );
a_c_err_olscb_PF = zeros( Nitr, Nk, 'mp' );
a_c_err_tlscb_PF = zeros( Nitr, Nk, 'mp' );
a_c_err_bezier_PF = zeros( Nitr, Nk, 'mp' );
a_c_err_lsqexp_PFX = zeros( Nitr, Nk, 'mp' );
a_c_err_olscb_PFX = zeros( Nitr, Nk, 'mp' );
a_c_err_tlscb_PFX = zeros( Nitr, Nk, 'mp' );
a_c_err_bezier_PFX = zeros( Nitr, Nk, 'mp' );
a_c_err_lsqexp_APX = zeros( Nitr, Nk, 'mp' );
a_c_err_olscb_APX = zeros( Nitr, Nk, 'mp' );
a_c_err_tlscb_APX = zeros( Nitr, Nk, 'mp' );
a_c_err_bezier_APX = zeros( Nitr, Nk, 'mp' );
a_c_err_lsqexp_APXF = zeros( Nitr, Nk, 'mp' );
a_c_err_olscb_APXF = zeros( Nitr, Nk, 'mp' );
a_c_err_tlscb_APXF = zeros( Nitr, Nk, 'mp' );
a_c_err_bezier_APXF = zeros( Nitr, Nk, 'mp' );
if ( b_do_extr )
    a_c_err_olscb_PF_extr = zeros( Nitr, Nk, 'mp' );
    a_c_err_tlscb_PF_extr = zeros( Nitr, Nk, 'mp' );
    a_c_err_olscb_PF_extr2 = zeros( Nitr, Nk, 'mp' );
    a_c_err_olscb_PFX_extr = zeros( Nitr, Nk, 'mp' );
    a_c_err_tlscb_PFX_extr = zeros( Nitr, Nk, 'mp' );
%    a_c_err_olscb_PFX_extr2 = zeros( Nitr, Nk, 'mp' );
    a_c_err_olscb_APX_extr = zeros( Nitr, Nk, 'mp' );
    a_c_err_tlscb_APX_extr = zeros( Nitr, Nk, 'mp' );
    a_c_err_olscb_APXF_extr = zeros( Nitr, Nk, 'mp' );
    a_c_err_tlscb_APXF_extr = zeros( Nitr, Nk, 'mp' );
end

v_c_nw_err_lsqexp_PF = zeros( Nitr, 1, 'mp' );
v_c_nw_err_olscb_PF = zeros( Nitr, 1, 'mp' );
v_c_nw_err_tlscb_PF = zeros( Nitr, 1, 'mp' );
v_c_nw_err_bezier_PF = zeros( Nitr, 1, 'mp' );
v_c_nw_err_lsqexp_PFX = zeros( Nitr, 1, 'mp' );
v_c_nw_err_olscb_PFX = zeros( Nitr, 1, 'mp' );
v_c_nw_err_tlscb_PFX = zeros( Nitr, 1, 'mp' );
v_c_nw_err_bezier_PFX = zeros( Nitr, 1, 'mp' );
v_c_nw_err_lsqexp_APX = zeros( Nitr, 1, 'mp' );
v_c_nw_err_olscb_APX = zeros( Nitr, 1, 'mp' );
v_c_nw_err_tlscb_APX = zeros( Nitr, 1, 'mp' );
v_c_nw_err_bezier_APX = zeros( Nitr, 1, 'mp' );
v_c_nw_err_lsqexp_APXF = zeros( Nitr, 1, 'mp' );
v_c_nw_err_olscb_APXF = zeros( Nitr, 1, 'mp' );
v_c_nw_err_tlscb_APXF = zeros( Nitr, 1, 'mp' );
v_c_nw_err_bezier_APXF = zeros( Nitr, 1, 'mp' );
if ( b_do_extr )
    v_c_nw_err_olscb_PF_extr = zeros( Nitr, 1, 'mp' );
    v_c_nw_err_tlscb_PF_extr = zeros( Nitr, 1, 'mp' );
    v_c_nw_err_olscb_PF_extr2 = zeros( Nitr, 1, 'mp' );
    v_c_nw_err_olscb_PFX_extr = zeros( Nitr, 1, 'mp' );
    v_c_nw_err_tlscb_PFX_extr = zeros( Nitr, 1, 'mp' );
%    v_c_nw_err_olscb_PFX_extr2 = zeros( Nitr, 1, 'mp' );
    v_c_nw_err_olscb_APX_extr = zeros( Nitr, 1, 'mp' );
    v_c_nw_err_tlscb_APX_extr = zeros( Nitr, 1, 'mp' );
    v_c_nw_err_olscb_APXF_extr = zeros( Nitr, 1, 'mp' );
    v_c_nw_err_tlscb_APXF_extr = zeros( Nitr, 1, 'mp' );
end




% Setup to find maximum and minimum v_c
v_c_lsqexp_PF_min = 1e50 * ones( 1, Nk );
v_c_lsqexp_PF_max = -1e50 * ones( 1, Nk );
v_c_olscb_PF_min = 1e50 * ones( 1, Nk );
v_c_olscb_PF_max = -1e50 * ones( 1, Nk );
v_c_tlscb_PF_min = 1e50 * ones( 1, Nk );
v_c_tlscb_PF_max = -1e50 * ones( 1, Nk );
v_c_bezier_PF_min = 1e50 * ones( 1, Nk );
v_c_bezier_PF_max = -1e50 * ones( 1, Nk );
v_c_lsqexp_PFX_min = 1e50 * ones( 1, Nk );
v_c_lsqexp_PFX_max = -1e50 * ones( 1, Nk );
v_c_olscb_PFX_min = 1e50 * ones( 1, Nk );
v_c_olscb_PFX_max = -1e50 * ones( 1, Nk );
v_c_tlscb_PFX_min = 1e50 * ones( 1, Nk );
v_c_tlscb_PFX_max = -1e50 * ones( 1, Nk );
v_c_bezier_PFX_min = 1e50 * ones( 1, Nk );
v_c_bezier_PFX_max = -1e50 * ones( 1, Nk );
v_c_lsqexp_APX_min = 1e50 * ones( 1, Nk );
v_c_lsqexp_APX_max = -1e50 * ones( 1, Nk );
v_c_olscb_APX_min = 1e50 * ones( 1, Nk );
v_c_olscb_APX_max = -1e50 * ones( 1, Nk );
v_c_tlscb_APX_min = 1e50 * ones( 1, Nk );
v_c_tlscb_APX_max = -1e50 * ones( 1, Nk );
v_c_bezier_APX_min = 1e50 * ones( 1, Nk );
v_c_bezier_APX_max = -1e50 * ones( 1, Nk );
v_c_lsqexp_APXF_min = 1e50 * ones( 1, Nk );
v_c_lsqexp_APXF_max = -1e50 * ones( 1, Nk );
v_c_olscb_APXF_min = 1e50 * ones( 1, Nk );
v_c_olscb_APXF_max = -1e50 * ones( 1, Nk );
v_c_tlscb_APXF_min = 1e50 * ones( 1, Nk );
v_c_tlscb_APXF_max = -1e50 * ones( 1, Nk );
v_c_bezier_APXF_min = 1e50 * ones( 1, Nk );
v_c_bezier_APXF_max = -1e50 * ones( 1, Nk );
if ( b_do_extr )
    v_c_olscb_PF_extr_min = 1e50 * ones( 1, Nk );
    v_c_olscb_PF_extr_max = -1e50 * ones( 1, Nk );
    v_c_tlscb_PF_extr_min = 1e50 * ones( 1, Nk );
    v_c_tlscb_PF_extr_max = -1e50 * ones( 1, Nk );
    v_c_olscb_PF_extr2_min = 1e50 * ones( 1, Nk );
    v_c_olscb_PF_extr2_max = -1e50 * ones( 1, Nk );
    v_c_olscb_PFX_extr_min = 1e50 * ones( 1, Nk );
    v_c_olscb_PFX_extr_max = -1e50 * ones( 1, Nk );
    v_c_tlscb_PFX_extr_min = 1e50 * ones( 1, Nk );
    v_c_tlscb_PFX_extr_max = -1e50 * ones( 1, Nk );
%     v_c_olscb_PFX_extr2_min = 1e50 * ones( 1, Nk );
%     v_c_olscb_PFX_extr2_max = -1e50 * ones( 1, Nk );
    v_c_olscb_APX_extr_min = 1e50 * ones( 1, Nk );
    v_c_olscb_APX_extr_max = -1e50 * ones( 1, Nk );
    v_c_tlscb_APX_extr_min = 1e50 * ones( 1, Nk );
    v_c_tlscb_APX_extr_max = -1e50 * ones( 1, Nk );
    v_c_olscb_APX_extr2_min = 1e50 * ones( 1, Nk );
    v_c_olscb_APX_extr2_max = -1e50 * ones( 1, Nk );
    v_c_olscb_APXF_extr_min = 1e50 * ones( 1, Nk );
    v_c_olscb_APXF_extr_max = -1e50 * ones( 1, Nk );
    v_c_tlscb_APXF_extr_min = 1e50 * ones( 1, Nk );
    v_c_tlscb_APXF_extr_max = -1e50 * ones( 1, Nk );
    v_c_olscb_APXF_extr2_min = 1e50 * ones( 1, Nk );
    v_c_olscb_APXF_extr2_max = -1e50 * ones( 1, Nk );
end


% ARGGH! Need this because of parfor loop
a_c_lsqexp_PF = zeros( Nitr, Nk );
a_c_olscb_PF = zeros( Nitr, Nk );
a_c_tlscb_PF = zeros( Nitr, Nk );
a_c_bezier_PF = zeros( Nitr, Nk );
a_c_lsqexp_PFX = zeros( Nitr, Nk );
a_c_olscb_PFX = zeros( Nitr, Nk );
a_c_tlscb_PFX = zeros( Nitr, Nk );
a_c_bezier_PFX = zeros( Nitr, Nk );
a_c_lsqexp_APX = zeros( Nitr, Nk );
a_c_olscb_APX = zeros( Nitr, Nk );
a_c_tlscb_APX = zeros( Nitr, Nk );
a_c_bezier_APX = zeros( Nitr, Nk );
a_c_lsqexp_APXF = zeros( Nitr, Nk );
a_c_olscb_APXF = zeros( Nitr, Nk );
a_c_tlscb_APXF = zeros( Nitr, Nk );
a_c_bezier_APXF = zeros( Nitr, Nk );
if ( b_do_extr )
    a_c_olscb_PF_extr = zeros( Nitr, Nk );
    a_c_tlscb_PF_extr = zeros( Nitr, Nk );
    a_c_olscb_PF_extr2 = zeros( Nitr, Nk );
    a_c_olscb_PFX_extr = zeros( Nitr, Nk );
    a_c_tlscb_PFX_extr = zeros( Nitr, Nk );
%    a_c_olscb_PFX_extr2 = zeros( Nitr, Nk );
    a_c_olscb_APX_extr = zeros( Nitr, Nk );
    a_c_tlscb_APX_extr = zeros( Nitr, Nk );
    a_c_olscb_APXF_extr = zeros( Nitr, Nk );
    a_c_tlscb_APXF_extr = zeros( Nitr, Nk );
end


ca_backup = cell( Nitr, 1 );

%ticBytes(gcp);
%warning( 'Not using parfor!' );
%for lp_itr=1:Nitr
parfor lp_itr=1:Nitr
    
    fprintf( '%d ', lp_itr );
    if ( 0 == mod( lp_itr, 30 ) ), fprintf( '\n' ); end
    
    % Do the calculation
    [ st_prtb_single, st_err ] = fn_ww__sim_shear_prtb__do_single( st_Dn, st_res, v_zs, v_zs_exact_U, v_zs_err_width, v_k, b_do_extr, st_p );
    ca_backup{lp_itr} = st_prtb_single;
    
    % Save for use in max min finding
    a_c_lsqexp_PF( lp_itr, : ) = st_prtb_single.v_c_lsqexp_PF;
    a_c_olscb_PF( lp_itr, : ) = st_prtb_single.v_c_olscb_PF;
    a_c_tlscb_PF( lp_itr, : ) = st_prtb_single.v_c_tlscb_PF;
    a_c_bezier_PF( lp_itr, : ) = st_prtb_single.v_c_bezier_PF;
    a_c_lsqexp_APX( lp_itr, : ) = st_prtb_single.v_c_lsqexp_APX;
    a_c_olscb_APX( lp_itr, : ) = st_prtb_single.v_c_olscb_APX;
    a_c_tlscb_APX( lp_itr, : ) = st_prtb_single.v_c_tlscb_APX;
    a_c_bezier_APX( lp_itr, : ) = st_prtb_single.v_c_bezier_APX;
    if ( b_do_extr )
        a_c_olscb_PF_extr( lp_itr, : ) = st_prtb_single.v_c_olscb_PF_extr;
        a_c_tlscb_PF_extr( lp_itr, : ) = st_prtb_single.v_c_tlscb_PF_extr;
        a_c_olscb_PF_extr2( lp_itr, : ) = st_prtb_single.v_c_olscb_PF_extr2;
        a_c_olscb_APX_extr( lp_itr, : ) = st_prtb_single.v_c_olscb_APX_extr;
        a_c_tlscb_APX_extr( lp_itr, : ) = st_prtb_single.v_c_tlscb_APX_extr;
    end    

    % Calculate componentwise relative errors
    a_c_err_lsqexp_PF( lp_itr, : ) = ( st_prtb_single.v_c_lsqexp_PF - v_c_ref ) ./ v_c_ref;    
    a_c_err_olscb_PF( lp_itr, : ) = ( st_prtb_single.v_c_olscb_PF - v_c_ref ) ./ v_c_ref;
    a_c_err_tlscb_PF( lp_itr, : ) = ( st_prtb_single.v_c_tlscb_PF - v_c_ref ) ./ v_c_ref;
    a_c_err_bezier_PF( lp_itr, : ) = ( st_prtb_single.v_c_bezier_PF - v_c_ref ) ./ v_c_ref;
    a_c_err_lsqexp_APX( lp_itr, : ) = ( st_prtb_single.v_c_lsqexp_APX - v_apx_ref ) ./ v_apx_ref;
    a_c_err_olscb_APX( lp_itr, : ) = ( st_prtb_single.v_c_olscb_APX - v_apx_ref ) ./ v_apx_ref;
    a_c_err_tlscb_APX( lp_itr, : ) = ( st_prtb_single.v_c_tlscb_APX - v_apx_ref ) ./ v_apx_ref;
    a_c_err_bezier_APX( lp_itr, : ) = ( st_prtb_single.v_c_bezier_APX - v_apx_ref ) ./ v_apx_ref;
    if ( b_do_extr )
        a_c_err_olscb_PF_extr( lp_itr, : ) = ( st_prtb_single.v_c_olscb_PF_extr - v_c_ref ) ./ v_c_ref;            
        a_c_err_tlscb_PF_extr( lp_itr, : ) = ( st_prtb_single.v_c_tlscb_PF_extr - v_c_ref ) ./ v_c_ref;
        a_c_err_olscb_PF_extr2( lp_itr, : ) = ( st_prtb_single.v_c_olscb_PF_extr2 - v_c_ref ) ./ v_c_ref;
        a_c_err_olscb_APX_extr( lp_itr, : ) = ( st_prtb_single.v_c_olscb_APX_extr - v_apx_ref ) ./ v_apx_ref;
        a_c_err_tlscb_APX_extr( lp_itr, : ) = ( st_prtb_single.v_c_tlscb_APX_extr - v_apx_ref ) ./ v_apx_ref;
    end

    % Calculate normwise relative errors
    norm_v_c_ref = norm( v_c_ref, inf );
    norm_v_apx_ref = norm( v_apx_ref, inf );
    v_c_nw_err_lsqexp_PF( lp_itr ) = norm( st_prtb_single.v_c_lsqexp_PF - v_c_ref, inf ) / norm_v_c_ref;
    v_c_nw_err_olscb_PF( lp_itr ) = norm( st_prtb_single.v_c_olscb_PF - v_c_ref, inf ) / norm_v_c_ref;
    v_c_nw_err_tlscb_PF( lp_itr ) = norm( st_prtb_single.v_c_tlscb_PF - v_c_ref, inf ) / norm_v_c_ref;
    v_c_nw_err_bezier_PF( lp_itr ) = norm( st_prtb_single.v_c_bezier_PF - v_c_ref, inf ) / norm_v_c_ref;
    v_c_nw_err_lsqexp_APX( lp_itr ) = norm( st_prtb_single.v_c_lsqexp_APX - v_apx_ref, inf ) / norm_v_apx_ref;
    v_c_nw_err_olscb_APX( lp_itr ) = norm( st_prtb_single.v_c_olscb_APX - v_apx_ref, inf ) / norm_v_apx_ref;
    v_c_nw_err_tlscb_APX( lp_itr ) = norm( st_prtb_single.v_c_tlscb_APX - v_apx_ref, inf ) / norm_v_apx_ref;
    v_c_nw_err_bezier_APX( lp_itr ) = norm( st_prtb_single.v_c_bezier_APX - v_apx_ref, inf ) / norm_v_apx_ref;
    if ( b_do_extr )
        v_c_nw_err_olscb_extr_PF( lp_itr ) = norm( st_prtb_single.v_c_olscb_PF_extr - v_c_ref, inf ) / norm_v_c_ref;
        v_c_nw_err_tlscb_extr_PF( lp_itr ) = norm( st_prtb_single.v_c_tlscb_PF_extr - v_c_ref, inf ) / norm_v_c_ref;
        v_c_nw_err_olscb_extr2_PF( lp_itr ) = norm( st_prtb_single.v_c_olscb_PF_extr2 - v_c_ref, inf ) / norm_v_c_ref;
        v_c_nw_err_olscb_extr_APX( lp_itr ) = norm( st_prtb_single.v_c_olscb_APX_extr - v_apx_ref, inf ) / norm_v_c_ref;
        v_c_nw_err_tlscb_extr_APX( lp_itr ) = norm( st_prtb_single.v_c_tlscb_APX_extr - v_apx_ref, inf ) / norm_v_c_ref;
    end 
        
    % Record error status
    if ( st_err.b_ok )                
        v_status( lp_itr ) = 1;        
    else       
        v_status( lp_itr ) = 0;        
    end    
    
end
fprintf( '\n\n' );
%tocBytes(gcp)





% Process any errors
v_pf_err_idxs = find( v_status == 0 )
if ( numel( v_pf_err_idxs > 0 ) )
    st_err.b_ok = false;
end
itr_done_pf = Nitr - numel( v_pf_err_idxs );
itr_done_apx = Nitr;




% First, create copy PF with crit layer entires replaced by APX
a_c_lsqexp_PFX = a_c_lsqexp_PF;
a_c_olscb_PFX = a_c_olscb_PF;
a_c_tlscb_PFX = a_c_tlscb_PF;
a_c_bezier_PFX = a_c_bezier_PF;
if ( b_do_extr )
    a_c_olscb_PFX_extr = a_c_olscb_PF_extr;
    a_c_tlscb_PFX_extr = a_c_tlscb_PF_extr;
%    a_c_olscb_PFX_extr2 = a_c_olscb_PF_extr2;
end
a_c_err_lsqexp_PFX = a_c_err_lsqexp_PF;
a_c_err_olscb_PFX = a_c_err_olscb_PF;
a_c_err_tlscb_PFX = a_c_err_tlscb_PF;
a_c_err_bezier_PFX = a_c_err_bezier_PF;
if ( b_do_extr )
    a_c_err_olscb_PFX_extr = a_c_err_olscb_PF_extr;
    a_c_err_tlscb_PFX_extr = a_c_err_tlscb_PF_extr;
%    a_c_err_olscb_PFX_extr2 = a_c_err_olscb_PF_extr2;
end
v_c_nw_err_lsqexp_PFX = v_c_nw_err_lsqexp_PF;
v_c_nw_err_olscb_PFX = v_c_nw_err_olscb_PF;
v_c_nw_err_tlscb_PFX = v_c_nw_err_tlscb_PF;
v_c_nw_err_bezier_PFX = v_c_nw_err_bezier_PF;
if ( b_do_extr )
    v_c_nw_err_olscb_extr_PFX = v_c_nw_err_olscb_extr_PF;
    v_c_nw_err_tlscb_extr_PFX = v_c_nw_err_tlscb_extr_PF;
    v_c_nw_err_olscb_extr2_PFX = v_c_nw_err_olscb_extr2_PF;
end
        
a_c_lsqexp_PFX( v_pf_err_idxs, : ) = a_c_lsqexp_APX( v_pf_err_idxs, : );
a_c_olscb_PFX( v_pf_err_idxs, : ) = a_c_olscb_APX( v_pf_err_idxs, : );
a_c_tlscb_PFX( v_pf_err_idxs, : ) = a_c_tlscb_APX( v_pf_err_idxs, : );
a_c_bezier_PFX( v_pf_err_idxs, : ) = a_c_bezier_APX( v_pf_err_idxs, : );
if ( b_do_extr )
    a_c_olscb_PFX_extr( v_pf_err_idxs, : ) = a_c_olscb_APX_extr( v_pf_err_idxs, : );
    a_c_tlscb_PFX_extr( v_pf_err_idxs, : ) = a_c_tlscb_APX_extr( v_pf_err_idxs, : );
%    a_c_olscb_PFX_extr2( v_pf_err_idxs, : ) = 0 * a_c_olscb_PFX_extr2( v_pf_err_idxs, : );
end
a_c_err_lsqexp_PFX( v_pf_err_idxs, : ) = a_c_err_lsqexp_APX( v_pf_err_idxs, : );
a_c_err_olscb_PFX( v_pf_err_idxs, : ) = a_c_err_olscb_APX( v_pf_err_idxs, : );
a_c_err_tlscb_PFX( v_pf_err_idxs, : ) = a_c_err_tlscb_APX( v_pf_err_idxs, : );
a_c_err_bezier_PFX( v_pf_err_idxs, : ) = a_c_err_bezier_APX( v_pf_err_idxs, : );
if ( b_do_extr )
    a_c_err_olscb_PFX_extr( v_pf_err_idxs, : ) = a_c_err_olscb_APX_extr( v_pf_err_idxs, : );
    a_c_err_tlscb_PFX_extr( v_pf_err_idxs, : ) = a_c_err_tlscb_APX_extr( v_pf_err_idxs, : );
%    a_c_err_olscb_PFX_extr2( v_pf_err_idxs, : ) = 0 * a_c_err_olscb_PFX_extr2( v_pf_err_idxs, : );
end
v_c_nw_err_lsqexp_PFX( v_pf_err_idxs ) = v_c_nw_err_lsqexp_APX( v_pf_err_idxs );
v_c_nw_err_olscb_PFX( v_pf_err_idxs ) = v_c_nw_err_olscb_APX( v_pf_err_idxs );
v_c_nw_err_tlscb_PFX( v_pf_err_idxs ) = v_c_nw_err_tlscb_APX( v_pf_err_idxs );
v_c_nw_err_bezier_PFX( v_pf_err_idxs ) = v_c_nw_err_bezier_APX( v_pf_err_idxs );
if ( b_do_extr )
    v_c_nw_err_olscb_extr_PFX( v_pf_err_idxs ) = v_c_nw_err_olscb_extr_APX( v_pf_err_idxs );
    v_c_nw_err_tlscb_extr_PFX( v_pf_err_idxs ) = v_c_nw_err_tlscb_extr_APX( v_pf_err_idxs );
%    v_c_nw_err_olscb_extr2_PFX( v_pf_err_idxs ) = 0 * v_c_nw_err_olscb_extr2_PFX( v_pf_err_idxs );
end




% Standard PF with crit layer entires removed entirely
a_c_lsqexp_PF( v_pf_err_idxs, : ) = [];
a_c_olscb_PF( v_pf_err_idxs, : ) = [];
a_c_tlscb_PF( v_pf_err_idxs, : ) = [];
a_c_bezier_PF( v_pf_err_idxs, : ) = [];
if ( b_do_extr )
    a_c_olscb_PF_extr( v_pf_err_idxs, : ) = [];
    a_c_tlscb_PF_extr( v_pf_err_idxs, : ) = [];
    a_c_olscb_PF_extr2( v_pf_err_idxs, : ) = [];
end
a_c_err_lsqexp_PF( v_pf_err_idxs, : ) = [];
a_c_err_olscb_PF( v_pf_err_idxs, : ) = [];
a_c_err_tlscb_PF( v_pf_err_idxs, : ) = [];
a_c_err_bezier_PF( v_pf_err_idxs, : ) = [];
if ( b_do_extr )
    a_c_err_olscb_PF_extr( v_pf_err_idxs, : ) = [];
    a_c_err_tlscb_PF_extr( v_pf_err_idxs, : ) = [];
    a_c_err_olscb_PF_extr2( v_pf_err_idxs, : ) = [];
end
v_c_nw_err_lsqexp_PF( v_pf_err_idxs ) = [];
v_c_nw_err_olscb_PF( v_pf_err_idxs ) = [];
v_c_nw_err_tlscb_PF( v_pf_err_idxs ) = [];
v_c_nw_err_bezier_PF( v_pf_err_idxs ) = [];
if ( b_do_extr )
    v_c_nw_err_olscb_extr_PF( v_pf_err_idxs ) = [];
    v_c_nw_err_tlscb_extr_PF( v_pf_err_idxs ) = [];
    v_c_nw_err_olscb_extr2_PF( v_pf_err_idxs ) = [];
end



% Create copy APX that is also filtered
a_c_lsqexp_APXF = a_c_lsqexp_APX;
a_c_olscb_APXF = a_c_olscb_APX;
a_c_tlscb_APXF = a_c_tlscb_APX;
a_c_bezier_APXF = a_c_bezier_APX;
if ( b_do_extr )
    a_c_olscb_APXF_extr = a_c_olscb_APX_extr;
    a_c_tlscb_APXF_extr = a_c_tlscb_APX_extr;
end    
a_c_err_lsqexp_APXF = a_c_err_lsqexp_APX;
a_c_err_olscb_APXF = a_c_err_olscb_APX;
a_c_err_tlscb_APXF = a_c_err_tlscb_APX;
a_c_err_bezier_APXF = a_c_err_bezier_APX;
if ( b_do_extr )
    a_c_err_olscb_APXF_extr = a_c_err_olscb_APX_extr;
    a_c_err_tlscb_APXF_extr = a_c_err_tlscb_APX_extr;
end
v_c_nw_err_lsqexp_APXF = v_c_nw_err_lsqexp_APX;
v_c_nw_err_olscb_APXF = v_c_nw_err_olscb_APX;
v_c_nw_err_tlscb_APXF = v_c_nw_err_tlscb_APX;
v_c_nw_err_bezier_APXF = v_c_nw_err_bezier_APX;
if ( b_do_extr )
    v_c_nw_err_olscb_extr_APXF = v_c_nw_err_olscb_extr_APX;
    v_c_nw_err_tlscb_extr_APXF = v_c_nw_err_tlscb_extr_APX;
end

a_c_lsqexp_APXF( v_pf_err_idxs, : ) = [];
a_c_olscb_APXF( v_pf_err_idxs, : ) = [];
a_c_tlscb_APXF( v_pf_err_idxs, : ) = [];
a_c_bezier_APXF( v_pf_err_idxs, : ) = [];
if ( b_do_extr )
    a_c_olscb_APXF_extr( v_pf_err_idxs, : ) = [];
    a_c_tlscb_APXF_extr( v_pf_err_idxs, : ) = [];
    %a_c_olscb_APXF_extr2( v_pf_err_idxs, : ) = [];
end
a_c_err_lsqexp_APXF( v_pf_err_idxs, : ) = [];
a_c_err_olscb_APXF( v_pf_err_idxs, : ) = [];
a_c_err_tlscb_APXF( v_pf_err_idxs, : ) = [];
a_c_err_bezier_APXF( v_pf_err_idxs, : ) = [];
if ( b_do_extr )
    a_c_err_olscb_APXF_extr( v_pf_err_idxs, : ) = [];
    a_c_err_tlscb_APXF_extr( v_pf_err_idxs, : ) = [];
    %a_c_err_olscb_APXF_extr2( v_pf_err_idxs, : ) = [];
end
v_c_nw_err_lsqexp_APXF( v_pf_err_idxs ) = [];
v_c_nw_err_olscb_APXF( v_pf_err_idxs ) = [];
v_c_nw_err_tlscb_APXF( v_pf_err_idxs ) = [];
v_c_nw_err_bezier_APXF( v_pf_err_idxs ) = [];
if ( b_do_extr )
    v_c_nw_err_olscb_extr_APXF( v_pf_err_idxs ) = [];
    v_c_nw_err_tlscb_extr_APXF( v_pf_err_idxs ) = [];
end








% Find min and max
v_c_lsqexp_PF_min = min( a_c_lsqexp_PF );
v_c_lsqexp_PF_max = max( a_c_lsqexp_PF );
v_c_olscb_PF_min = min( a_c_olscb_PF );
v_c_olscb_PF_max = max( a_c_olscb_PF );
v_c_tlscb_PF_min = min( a_c_tlscb_PF );
v_c_tlscb_PF_max = max( a_c_tlscb_PF );
v_c_bezier_PF_min = min( a_c_bezier_PF );
v_c_bezier_PF_max = max( a_c_bezier_PF );
v_c_lsqexp_PFX_min = min( a_c_lsqexp_PFX );
v_c_lsqexp_PFX_max = max( a_c_lsqexp_PFX );
v_c_olscb_PFX_min = min( a_c_olscb_PFX );
v_c_olscb_PFX_max = max( a_c_olscb_PFX );
v_c_tlscb_PFX_min = min( a_c_tlscb_PFX );
v_c_tlscb_PFX_max = max( a_c_tlscb_PFX );
v_c_bezier_PFX_min = min( a_c_bezier_PFX );
v_c_bezier_PFX_max = max( a_c_bezier_PFX );
v_c_lsqexp_APX_min = min( a_c_lsqexp_APX );
v_c_lsqexp_APX_max = max( a_c_lsqexp_APX );
v_c_olscb_APX_min = min( a_c_olscb_APX );
v_c_olscb_APX_max = max( a_c_olscb_APX );
v_c_tlscb_APX_min = min( a_c_tlscb_APX );
v_c_tlscb_APX_max = max( a_c_tlscb_APX );
v_c_bezier_APX_min = min( a_c_bezier_APX );
v_c_bezier_APX_max = max( a_c_bezier_APX );
v_c_lsqexp_APXF_min = min( a_c_lsqexp_APXF );
v_c_lsqexp_APXF_max = max( a_c_lsqexp_APXF );
v_c_olscb_APXF_min = min( a_c_olscb_APXF );
v_c_olscb_APXF_max = max( a_c_olscb_APXF );
v_c_tlscb_APXF_min = min( a_c_tlscb_APXF );
v_c_tlscb_APXF_max = max( a_c_tlscb_APXF );
v_c_bezier_APXF_min = min( a_c_bezier_APXF );
v_c_bezier_APXF_max = max( a_c_bezier_APXF );
if ( b_do_extr )
    v_c_olscb_PF_extr_min = min( a_c_olscb_PF_extr );
    v_c_olscb_PF_extr_max = max( a_c_olscb_PF_extr );
    v_c_tlscb_PF_extr_min = min( a_c_tlscb_PF_extr );
    v_c_tlscb_PF_extr_max = max( a_c_tlscb_PF_extr );
    v_c_olscb_PF_extr2_min = min( a_c_olscb_PF_extr2 );
    v_c_olscb_PF_extr2_max = max( a_c_olscb_PF_extr2 );
    v_c_olscb_PFX_extr_min = min( a_c_olscb_PFX_extr );
    v_c_olscb_PFX_extr_max = max( a_c_olscb_PFX_extr );
    v_c_tlscb_PFX_extr_min = min( a_c_tlscb_PFX_extr );
    v_c_tlscb_PFX_extr_max = max( a_c_tlscb_PFX_extr );
%     v_c_olscb_PFX_extr2_min = min( a_c_olscb_PFX_extr2 );
%     v_c_olscb_PFX_extr2_max = max( a_c_olscb_PFX_extr2 );        
    v_c_olscb_APX_extr_min = min( a_c_olscb_APX_extr );
    v_c_olscb_APX_extr_max = max( a_c_olscb_APX_extr );
    v_c_tlscb_APX_extr_min = min( a_c_tlscb_APX_extr );
    v_c_tlscb_APX_extr_max = max( a_c_tlscb_APX_extr );
    v_c_olscb_APXF_extr_min = min( a_c_olscb_APXF_extr );
    v_c_olscb_APXF_extr_max = max( a_c_olscb_APXF_extr );
    v_c_tlscb_APXF_extr_min = min( a_c_tlscb_APXF_extr );
    v_c_tlscb_APXF_extr_max = max( a_c_tlscb_APXF_extr );
end



% Calculate statistical moments, must use high precision arithmetic or such
% like when available
mean_err_lsqexp_PF = mean( a_c_err_lsqexp_PF );
mean_err_olscb_PF = mean( a_c_err_olscb_PF );
mean_err_tlscb_PF = mean( a_c_err_tlscb_PF );
mean_err_bezier_PF = mean( a_c_err_bezier_PF );
mean_err_lsqexp_PFX = mean( a_c_err_lsqexp_PFX );
mean_err_olscb_PFX = mean( a_c_err_olscb_PFX );
mean_err_tlscb_PFX = mean( a_c_err_tlscb_PFX );
mean_err_bezier_PFX = mean( a_c_err_bezier_PFX );
mean_err_lsqexp_APX = mean( a_c_err_lsqexp_APX );
mean_err_olscb_APX = mean( a_c_err_olscb_APX );
mean_err_tlscb_APX = mean( a_c_err_tlscb_APX );
mean_err_bezier_APX = mean( a_c_err_bezier_APX );
mean_err_lsqexp_APXF = mean( a_c_err_lsqexp_APXF );
mean_err_olscb_APXF = mean( a_c_err_olscb_APXF );
mean_err_tlscb_APXF = mean( a_c_err_tlscb_APXF );
mean_err_bezier_APXF = mean( a_c_err_bezier_APXF );

std_err_lsqexp_PF = std( a_c_err_lsqexp_PF );
std_err_olscb_PF = std( a_c_err_olscb_PF );
std_err_tlscb_PF = std( a_c_err_tlscb_PF );
std_err_bezier_PF = std( a_c_err_bezier_PF );
std_err_lsqexp_PFX = std( a_c_err_lsqexp_PFX );
std_err_olscb_PFX = std( a_c_err_olscb_PFX );
std_err_tlscb_PFX = std( a_c_err_tlscb_PFX );
std_err_bezier_PFX = std( a_c_err_bezier_PFX );
std_err_lsqexp_APX = std( a_c_err_lsqexp_APX );
std_err_olscb_APX = std( a_c_err_olscb_APX );
std_err_tlscb_APX = std( a_c_err_tlscb_APX );
std_err_bezier_APX = std( a_c_err_bezier_APX );
std_err_lsqexp_APXF = std( a_c_err_lsqexp_APXF );
std_err_olscb_APXF = std( a_c_err_olscb_APXF );
std_err_tlscb_APXF = std( a_c_err_tlscb_APXF );
std_err_bezier_APXF = std( a_c_err_bezier_APXF );

skew_err_lsqexp_PF = skewness( double( a_c_err_lsqexp_PF ) );
skew_err_olscb_PF = skewness( double( a_c_err_olscb_PF ) );
skew_err_tlscb_PF = skewness( double( a_c_err_tlscb_PF ) );
skew_err_bezier_PF = skewness( double( a_c_err_bezier_PF ) );
skew_err_lsqexp_PFX = skewness( double( a_c_err_lsqexp_PFX ) );
skew_err_olscb_PFX = skewness( double( a_c_err_olscb_PFX ) );
skew_err_tlscb_PFX = skewness( double( a_c_err_tlscb_PFX ) );
skew_err_bezier_PFX = skewness( double( a_c_err_bezier_PFX ) );
skew_err_lsqexp_APX = skewness( double( a_c_err_lsqexp_APX ) );
skew_err_olscb_APX = skewness( double( a_c_err_olscb_APX ) );
skew_err_tlscb_APX = skewness( double( a_c_err_tlscb_APX ) );
skew_err_bezier_APX = skewness( double( a_c_err_bezier_APX ) );
skew_err_lsqexp_APXF = skewness( double( a_c_err_lsqexp_APXF ) );
skew_err_olscb_APXF = skewness( double( a_c_err_olscb_APXF ) );
skew_err_tlscb_APXF = skewness( double( a_c_err_tlscb_APXF ) );
skew_err_bezier_APXF = skewness( double( a_c_err_bezier_APXF ) );

kurt_err_lsqexp_PF = kurtosis( double( a_c_err_lsqexp_PF ) ) - 3;
kurt_err_olscb_PF = kurtosis( double( a_c_err_olscb_PF ) ) - 3;
kurt_err_tlscb_PF = kurtosis( double( a_c_err_tlscb_PF ) ) - 3;
kurt_err_bezier_PF = kurtosis( double( a_c_err_bezier_PF ) ) - 3;
kurt_err_lsqexp_PFX = kurtosis( double( a_c_err_lsqexp_PFX ) ) - 3;
kurt_err_olscb_PFX = kurtosis( double( a_c_err_olscb_PFX ) ) - 3;
kurt_err_tlscb_PFX = kurtosis( double( a_c_err_tlscb_PFX ) ) - 3;
kurt_err_bezier_PFX = kurtosis( double( a_c_err_bezier_PFX ) ) - 3;
kurt_err_lsqexp_APX = kurtosis( double( a_c_err_lsqexp_APX ) ) - 3;
kurt_err_olscb_APX = kurtosis( double( a_c_err_olscb_APX ) ) - 3;
kurt_err_tlscb_APX = kurtosis( double( a_c_err_tlscb_APX ) ) - 3;
kurt_err_bezier_APX = kurtosis( double( a_c_err_bezier_APX ) ) - 3;
kurt_err_lsqexp_APXF = kurtosis( double( a_c_err_lsqexp_APXF ) ) - 3;
kurt_err_olscb_APXF = kurtosis( double( a_c_err_olscb_APXF ) ) - 3;
kurt_err_tlscb_APXF = kurtosis( double( a_c_err_tlscb_APXF ) ) - 3;
kurt_err_bezier_APXF = kurtosis( double( a_c_err_bezier_APXF ) ) - 3;
if ( b_do_extr )
    mean_err_olscb_PF_extr = mean( a_c_err_olscb_PF_extr );
    mean_err_tlscb_PF_extr = mean( a_c_err_tlscb_PF_extr );
    mean_err_olscb_PF_extr2 = mean( a_c_err_olscb_PF_extr2 );
    mean_err_olscb_PFX_extr = mean( a_c_err_olscb_PFX_extr );
    mean_err_tlscb_PFX_extr = mean( a_c_err_tlscb_PFX_extr );
%     mean_err_olscb_PFX_extr2 = mean( a_c_err_olscb_PFX_extr2 );
    mean_err_olscb_APX_extr = mean( a_c_err_olscb_APX_extr );
    mean_err_tlscb_APX_extr = mean( a_c_err_tlscb_APX_extr );
    mean_err_olscb_APXF_extr = mean( a_c_err_olscb_APXF_extr );
    mean_err_tlscb_APXF_extr = mean( a_c_err_tlscb_APXF_extr );
    
    std_err_olscb_PF_extr = std( a_c_err_olscb_PF_extr );
    std_err_tlscb_PF_extr = std( a_c_err_tlscb_PF_extr );
    std_err_olscb_PF_extr2 = std( a_c_err_olscb_PF_extr2 );
    std_err_olscb_PFX_extr = std( a_c_err_olscb_PFX_extr );
    std_err_tlscb_PFX_extr = std( a_c_err_tlscb_PFX_extr );
%     std_err_olscb_PFX_extr2 = std( a_c_err_olscb_PFX_extr2 );
    std_err_olscb_APX_extr = std( a_c_err_olscb_APX_extr );
    std_err_tlscb_APX_extr = std( a_c_err_tlscb_APX_extr );
    std_err_olscb_APXF_extr = std( a_c_err_olscb_APXF_extr );
    std_err_tlscb_APXF_extr = std( a_c_err_tlscb_APXF_extr );
    
    skew_err_olscb_PF_extr = skewness( double( a_c_err_olscb_PF_extr ) );
    skew_err_tlscb_PF_extr = skewness( double( a_c_err_tlscb_PF_extr ) );
    skew_err_olscb_PF_extr2 = skewness( double( a_c_err_olscb_PF_extr2 ) );
    skew_err_olscb_PFX_extr = skewness( double( a_c_err_olscb_PFX_extr ) );
    skew_err_tlscb_PFX_extr = skewness( double( a_c_err_tlscb_PFX_extr ) );
%     skew_err_olscb_PFX_extr2 = skewness( double( a_c_err_olscb_PFX_extr2 ) );
    skew_err_olscb_APX_extr = skewness( double( a_c_err_olscb_APX_extr ) );
    skew_err_tlscb_APX_extr = skewness( double( a_c_err_tlscb_APX_extr ) );
    skew_err_olscb_APXF_extr = skewness( double( a_c_err_olscb_APXF_extr ) );
    skew_err_tlscb_APXF_extr = skewness( double( a_c_err_tlscb_APXF_extr ) );
    
    kurt_err_olscb_PF_extr = kurtosis( double( a_c_err_olscb_PF_extr ) ) - 3;
    kurt_err_tlscb_PF_extr = kurtosis( double( a_c_err_tlscb_PF_extr ) ) - 3;
    kurt_err_olscb_PF_extr2 = kurtosis( double( a_c_err_olscb_PF_extr2 ) ) - 3;
    kurt_err_olscb_PFX_extr = kurtosis( double( a_c_err_olscb_PFX_extr ) ) - 3;
    kurt_err_tlscb_PFX_extr = kurtosis( double( a_c_err_tlscb_PFX_extr ) ) - 3;
%     kurt_err_olscb_PFX_extr2 = kurtosis( double( a_c_err_olscb_PFX_extr2 ) ) - 3;
    kurt_err_olscb_APX_extr = kurtosis( double( a_c_err_olscb_APX_extr ) ) - 3;
    kurt_err_tlscb_APX_extr = kurtosis( double( a_c_err_tlscb_APX_extr ) ) - 3;
    kurt_err_olscb_APXF_extr = kurtosis( double( a_c_err_olscb_APXF_extr ) ) - 3;
    kurt_err_tlscb_APXF_extr = kurtosis( double( a_c_err_tlscb_APXF_extr ) ) - 3;
end




% Calculate normwise stat quantities
mean_nw_err_lsqexp_PF = mean( v_c_nw_err_lsqexp_PF );
mean_nw_err_olscb_PF = mean( v_c_nw_err_olscb_PF );
mean_nw_err_tlscb_PF = mean( v_c_nw_err_tlscb_PF );
mean_nw_err_bezier_PF = mean( v_c_nw_err_bezier_PF );
mean_nw_err_lsqexp_PFX = mean( v_c_nw_err_lsqexp_PFX );
mean_nw_err_olscb_PFX = mean( v_c_nw_err_olscb_PFX );
mean_nw_err_tlscb_PFX = mean( v_c_nw_err_tlscb_PFX );
mean_nw_err_bezier_PFX = mean( v_c_nw_err_bezier_PFX );
mean_nw_err_lsqexp_APX = mean( v_c_nw_err_lsqexp_APX );
mean_nw_err_olscb_APX = mean( v_c_nw_err_olscb_APX );
mean_nw_err_tlscb_APX = mean( v_c_nw_err_tlscb_APX );
mean_nw_err_bezier_APX = mean( v_c_nw_err_bezier_APX );
mean_nw_err_lsqexp_APXF = mean( v_c_nw_err_lsqexp_APXF );
mean_nw_err_olscb_APXF = mean( v_c_nw_err_olscb_APXF );
mean_nw_err_tlscb_APXF = mean( v_c_nw_err_tlscb_APXF );
mean_nw_err_bezier_APXF = mean( v_c_nw_err_bezier_APXF );

std_nw_err_lsqexp_PF = std( v_c_nw_err_lsqexp_PF );
std_nw_err_olscb_PF = std( v_c_nw_err_olscb_PF );
std_nw_err_tlscb_PF = std( v_c_nw_err_tlscb_PF );
std_nw_err_bezier_PF = std( v_c_nw_err_bezier_PF );
std_nw_err_lsqexp_PFX = std( v_c_nw_err_lsqexp_PFX );
std_nw_err_olscb_PFX = std( v_c_nw_err_olscb_PFX );
std_nw_err_tlscb_PFX = std( v_c_nw_err_tlscb_PFX );
std_nw_err_bezier_PFX = std( v_c_nw_err_bezier_PFX );
std_nw_err_lsqexp_APX = std( v_c_nw_err_lsqexp_APX );
std_nw_err_olscb_APX = std( v_c_nw_err_olscb_APX );
std_nw_err_tlscb_APX = std( v_c_nw_err_tlscb_APX );
std_nw_err_bezier_APX = std( v_c_nw_err_bezier_APX );
std_nw_err_lsqexp_APXF = std( v_c_nw_err_lsqexp_APXF );
std_nw_err_olscb_APXF = std( v_c_nw_err_olscb_APXF );
std_nw_err_tlscb_APXF = std( v_c_nw_err_tlscb_APXF );
std_nw_err_bezier_APXF = std( v_c_nw_err_bezier_APXF );

skew_nw_err_lsqexp_PF = skewness( double( v_c_nw_err_lsqexp_PF ) );
skew_nw_err_olscb_PF = skewness( double( v_c_nw_err_olscb_PF ) );
skew_nw_err_tlscb_PF = skewness( double( v_c_nw_err_tlscb_PF ) );
skew_nw_err_bezier_PF = skewness( double( v_c_nw_err_bezier_PF ) );
skew_nw_err_lsqexp_PFX = skewness( double( v_c_nw_err_lsqexp_PFX ) );
skew_nw_err_olscb_PFX = skewness( double( v_c_nw_err_olscb_PFX ) );
skew_nw_err_tlscb_PFX = skewness( double( v_c_nw_err_tlscb_PFX ) );
skew_nw_err_bezier_PFX = skewness( double( v_c_nw_err_bezier_PFX ) );
skew_nw_err_lsqexp_APX = skewness( double( v_c_nw_err_lsqexp_APX ) );
skew_nw_err_olscb_APX = skewness( double( v_c_nw_err_olscb_APX ) );
skew_nw_err_tlscb_APX = skewness( double( v_c_nw_err_tlscb_APX ) );
skew_nw_err_bezier_APX = skewness( double( v_c_nw_err_bezier_APX ) );
skew_nw_err_lsqexp_APXF = skewness( double( v_c_nw_err_lsqexp_APXF ) );
skew_nw_err_olscb_APXF = skewness( double( v_c_nw_err_olscb_APXF ) );
skew_nw_err_tlscb_APXF = skewness( double( v_c_nw_err_tlscb_APXF ) );
skew_nw_err_bezier_APXF = skewness( double( v_c_nw_err_bezier_APXF ) );

kurt_nw_err_lsqexp_PF = kurtosis( double( v_c_nw_err_lsqexp_PF ) );
kurt_nw_err_olscb_PF = kurtosis( double( v_c_nw_err_olscb_PF ) );
kurt_nw_err_tlscb_PF = kurtosis( double( v_c_nw_err_tlscb_PF ) );
kurt_nw_err_bezier_PF = kurtosis( double( v_c_nw_err_bezier_PF ) );
kurt_nw_err_lsqexp_PFX = kurtosis( double( v_c_nw_err_lsqexp_PFX ) );
kurt_nw_err_olscb_PFX = kurtosis( double( v_c_nw_err_olscb_PFX ) );
kurt_nw_err_tlscb_PFX = kurtosis( double( v_c_nw_err_tlscb_PFX ) );
kurt_nw_err_bezier_PFX = kurtosis( double( v_c_nw_err_bezier_PFX ) );
kurt_nw_err_lsqexp_APX = kurtosis( double( v_c_nw_err_lsqexp_APX ) );
kurt_nw_err_olscb_APX = kurtosis( double( v_c_nw_err_olscb_APX ) );
kurt_nw_err_tlscb_APX = kurtosis( double( v_c_nw_err_tlscb_APX ) );
kurt_nw_err_bezier_APX = kurtosis( double( v_c_nw_err_bezier_APX ) );
kurt_nw_err_lsqexp_APXF = kurtosis( double( v_c_nw_err_lsqexp_APXF ) );
kurt_nw_err_olscb_APXF = kurtosis( double( v_c_nw_err_olscb_APXF ) );
kurt_nw_err_tlscb_APXF = kurtosis( double( v_c_nw_err_tlscb_APXF ) );
kurt_nw_err_bezier_APXF = kurtosis( double( v_c_nw_err_bezier_APXF ) );

if ( b_do_extr )
    mean_nw_err_olscb_PF_extr = mean( v_c_nw_err_olscb_PF_extr );
    mean_nw_err_tlscb_PF_extr = mean( v_c_nw_err_tlscb_PF_extr );
    mean_nw_err_olscb_PF_extr2 = mean( v_c_nw_err_olscb_PF_extr2 );
    mean_nw_err_olscb_PFX_extr = mean( v_c_nw_err_olscb_PFX_extr );
    mean_nw_err_tlscb_PFX_extr = mean( v_c_nw_err_tlscb_PFX_extr );
%     mean_nw_err_olscb_PFX_extr2 = mean( v_c_nw_err_olscb_PFX_extr2 );
    mean_nw_err_olscb_APX_extr = mean( v_c_nw_err_olscb_APX_extr );
    mean_nw_err_tlscb_APX_extr = mean( v_c_nw_err_tlscb_APX_extr );
    mean_nw_err_olscb_APXF_extr = mean( v_c_nw_err_olscb_APXF_extr );
    mean_nw_err_tlscb_APXF_extr = mean( v_c_nw_err_tlscb_APXF_extr );
    
    std_nw_err_olscb_PF_extr = std( v_c_nw_err_olscb_PF_extr );
    std_nw_err_tlscb_PF_extr = std( v_c_nw_err_tlscb_PF_extr );
    std_nw_err_olscb_PF_extr2 = std( v_c_nw_err_olscb_PF_extr2 );
    std_nw_err_olscb_PFX_extr = std( v_c_nw_err_olscb_PFX_extr );
    std_nw_err_tlscb_PFX_extr = std( v_c_nw_err_tlscb_PFX_extr );
%     std_nw_err_olscb_PFX_extr2 = std( v_c_nw_err_olscb_PFX_extr2 );
    std_nw_err_olscb_APX_extr = std( v_c_nw_err_olscb_APX_extr );
    std_nw_err_tlscb_APX_extr = std( v_c_nw_err_tlscb_APX_extr );
    std_nw_err_olscb_APXF_extr = std( v_c_nw_err_olscb_APXF_extr );
    std_nw_err_tlscb_APXF_extr = std( v_c_nw_err_tlscb_APXF_extr );
    
    skew_nw_err_olscb_PF_extr = skewness( double( v_c_nw_err_olscb_PF_extr ) );
    skew_nw_err_tlscb_PF_extr = skewness( double( v_c_nw_err_tlscb_PF_extr ) );
    skew_nw_err_olscb_PF_extr2 = skewness( double( v_c_nw_err_olscb_PF_extr2 ) );
    skew_nw_err_olscb_PFX_extr = skewness( double( v_c_nw_err_olscb_PFX_extr ) );
    skew_nw_err_tlscb_PFX_extr = skewness( double( v_c_nw_err_tlscb_PFX_extr ) );
%    skew_nw_err_olscb_PFX_extr2 = skewness( double( v_c_nw_err_olscb_PFX_extr2 ) );
    skew_nw_err_olscb_APX_extr = skewness( double( v_c_nw_err_olscb_APX_extr ) );
    skew_nw_err_tlscb_APX_extr = skewness( double( v_c_nw_err_tlscb_APX_extr ) );
    skew_nw_err_olscb_APXF_extr = skewness( double( v_c_nw_err_olscb_APXF_extr ) );
    skew_nw_err_tlscb_APXF_extr = skewness( double( v_c_nw_err_tlscb_APXF_extr ) );
    
    kurt_nw_err_olscb_PF_extr = kurtosis( double( v_c_nw_err_olscb_PF_extr ) );
    kurt_nw_err_tlscb_PF_extr = kurtosis( double( v_c_nw_err_tlscb_PF_extr ) );
    kurt_nw_err_olscb_PF_extr2 = kurtosis( double( v_c_nw_err_olscb_PF_extr2 ) );
    kurt_nw_err_olscb_PFX_extr = kurtosis( double( v_c_nw_err_olscb_PFX_extr ) );
    kurt_nw_err_tlscb_PFX_extr = kurtosis( double( v_c_nw_err_tlscb_PFX_extr ) );
%     kurt_nw_err_olscb_PFX_extr2 = kurtosis( double( v_c_nw_err_olscb_PFX_extr2 ) );
    kurt_nw_err_olscb_APX_extr = kurtosis( double( v_c_nw_err_olscb_APX_extr ) );
    kurt_nw_err_tlscb_APX_extr = kurtosis( double( v_c_nw_err_tlscb_APX_extr ) );
    kurt_nw_err_olscb_APXF_extr = kurtosis( double( v_c_nw_err_olscb_APXF_extr ) );
    kurt_nw_err_tlscb_APXF_extr = kurtosis( double( v_c_nw_err_tlscb_APXF_extr ) );
end


% Calculate percentile quantities
a_prctl_mean_err_lsqexp_PF = prctile( double( a_c_err_lsqexp_PF ), v_percentile );
a_prctl_mean_err_olscb_PF = prctile( double( a_c_err_olscb_PF ), v_percentile );
a_prctl_mean_err_tlscb_PF = prctile( double( a_c_err_tlscb_PF ), v_percentile );
a_prctl_mean_err_bezier_PF = prctile( double( a_c_err_bezier_PF ), v_percentile );
a_prctl_mean_err_lsqexp_PFX = prctile( double( a_c_err_lsqexp_PFX ), v_percentile );
a_prctl_mean_err_olscb_PFX = prctile( double( a_c_err_olscb_PFX ), v_percentile );
a_prctl_mean_err_tlscb_PFX = prctile( double( a_c_err_tlscb_PFX ), v_percentile );
a_prctl_mean_err_bezier_PFX = prctile( double( a_c_err_bezier_PFX ), v_percentile );
a_prctl_mean_err_lsqexp_APX = prctile( double( a_c_err_lsqexp_APX ), v_percentile );
a_prctl_mean_err_olscb_APX = prctile( double( a_c_err_olscb_APX ), v_percentile );
a_prctl_mean_err_tlscb_APX = prctile( double( a_c_err_tlscb_APX ), v_percentile );
a_prctl_mean_err_bezier_APX = prctile( double( a_c_err_bezier_APX ), v_percentile );
a_prctl_mean_err_lsqexp_APXF = prctile( double( a_c_err_lsqexp_APXF ), v_percentile );
a_prctl_mean_err_olscb_APXF = prctile( double( a_c_err_olscb_APXF ), v_percentile );
a_prctl_mean_err_tlscb_APXF = prctile( double( a_c_err_tlscb_APXF ), v_percentile );
a_prctl_mean_err_bezier_APXF = prctile( double( a_c_err_bezier_APXF ), v_percentile );
if ( b_do_extr )
    a_prctl_mean_err_olscb_PF_extr = prctile( double( a_c_err_olscb_PF_extr ), v_percentile );
    a_prctl_mean_err_tlscb_PF_extr = prctile( double( a_c_err_tlscb_PF_extr ), v_percentile );
    a_prctl_mean_err_olscb_PF_extr2 = prctile( double( a_c_err_olscb_PF_extr2 ), v_percentile );
    a_prctl_mean_err_olscb_PFX_extr = prctile( double( a_c_err_olscb_PFX_extr ), v_percentile );
    a_prctl_mean_err_tlscb_PFX_extr = prctile( double( a_c_err_tlscb_PFX_extr ), v_percentile );
%     a_prctl_mean_err_olscb_PFX_extr2 = prctile( double( a_c_err_olscb_PFX_extr2 ), v_percentile );
    a_prctl_mean_err_olscb_APX_extr = prctile( double( a_c_err_olscb_APX_extr ), v_percentile );
    a_prctl_mean_err_tlscb_APX_extr = prctile( double( a_c_err_tlscb_APX_extr ), v_percentile );
    a_prctl_mean_err_olscb_APXF_extr = prctile( double( a_c_err_olscb_APXF_extr ), v_percentile );
    a_prctl_mean_err_tlscb_APXF_extr = prctile( double( a_c_err_tlscb_APXF_extr ), v_percentile );
end


% Dump into struct
st_prtb_results = struct;

st_prtb_results.mean_err_lsqexp_PF = mean_err_lsqexp_PF;
st_prtb_results.mean_err_olscb_PF = mean_err_olscb_PF;
st_prtb_results.mean_err_tlscb_PF = mean_err_tlscb_PF;
st_prtb_results.mean_err_bezier_PF = mean_err_bezier_PF;
st_prtb_results.mean_err_lsqexp_PFX = mean_err_lsqexp_PFX;
st_prtb_results.mean_err_olscb_PFX = mean_err_olscb_PFX;
st_prtb_results.mean_err_tlscb_PFX = mean_err_tlscb_PFX;
st_prtb_results.mean_err_bezier_PFX = mean_err_bezier_PFX;
st_prtb_results.mean_err_lsqexp_APX = mean_err_lsqexp_APX;
st_prtb_results.mean_err_olscb_APX = mean_err_olscb_APX;
st_prtb_results.mean_err_tlscb_APX = mean_err_tlscb_APX;
st_prtb_results.mean_err_bezier_APX = mean_err_bezier_APX;
st_prtb_results.mean_err_lsqexp_APXF = mean_err_lsqexp_APXF;
st_prtb_results.mean_err_olscb_APXF = mean_err_olscb_APXF;
st_prtb_results.mean_err_tlscb_APXF = mean_err_tlscb_APXF;
st_prtb_results.mean_err_bezier_APXF = mean_err_bezier_APXF;

st_prtb_results.std_err_lsqexp_PF = std_err_lsqexp_PF;
st_prtb_results.std_err_olscb_PF = std_err_olscb_PF;
st_prtb_results.std_err_tlscb_PF = std_err_tlscb_PF;
st_prtb_results.std_err_bezier_PF = std_err_bezier_PF;
st_prtb_results.std_err_lsqexp_PFX = std_err_lsqexp_PFX;
st_prtb_results.std_err_olscb_PFX = std_err_olscb_PFX;
st_prtb_results.std_err_tlscb_PFX = std_err_tlscb_PFX;
st_prtb_results.std_err_bezier_PFX = std_err_bezier_PFX;
st_prtb_results.std_err_lsqexp_APX = std_err_lsqexp_APX;
st_prtb_results.std_err_olscb_APX = std_err_olscb_APX;
st_prtb_results.std_err_tlscb_APX = std_err_tlscb_APX;
st_prtb_results.std_err_bezier_APX = std_err_bezier_APX;
st_prtb_results.std_err_lsqexp_APXF = std_err_lsqexp_APXF;
st_prtb_results.std_err_olscb_APXF = std_err_olscb_APXF;
st_prtb_results.std_err_tlscb_APXF = std_err_tlscb_APXF;
st_prtb_results.std_err_bezier_APXF = std_err_bezier_APXF;

st_prtb_results.skew_err_lsqexp_PF = skew_err_lsqexp_PF;
st_prtb_results.skew_err_olscb_PF = skew_err_olscb_PF;
st_prtb_results.skew_err_tlscb_PF = skew_err_tlscb_PF;
st_prtb_results.skew_err_bezier_PF = skew_err_bezier_PF;
st_prtb_results.skew_err_lsqexp_PFX = skew_err_lsqexp_PFX;
st_prtb_results.skew_err_olscb_PFX = skew_err_olscb_PFX;
st_prtb_results.skew_err_tlscb_PFX = skew_err_tlscb_PFX;
st_prtb_results.skew_err_bezier_PFX = skew_err_bezier_PFX;
st_prtb_results.skew_err_lsqexp_APX = skew_err_lsqexp_APX;
st_prtb_results.skew_err_olscb_APX = skew_err_olscb_APX;
st_prtb_results.skew_err_tlscb_APX = skew_err_tlscb_APX;
st_prtb_results.skew_err_bezier_APX = skew_err_bezier_APX;
st_prtb_results.skew_err_lsqexp_APXF = skew_err_lsqexp_APXF;
st_prtb_results.skew_err_olscb_APXF = skew_err_olscb_APXF;
st_prtb_results.skew_err_tlscb_APXF = skew_err_tlscb_APXF;
st_prtb_results.skew_err_bezier_APXF = skew_err_bezier_APXF;

st_prtb_results.kurt_err_lsqexp_PF = kurt_err_lsqexp_PF;
st_prtb_results.kurt_err_olscb_PF = kurt_err_olscb_PF;
st_prtb_results.kurt_err_tlscb_PF = kurt_err_tlscb_PF;
st_prtb_results.kurt_err_bezier_PF = kurt_err_bezier_PF;
st_prtb_results.kurt_err_lsqexp_PFX = kurt_err_lsqexp_PFX;
st_prtb_results.kurt_err_olscb_PFX = kurt_err_olscb_PFX;
st_prtb_results.kurt_err_tlscb_PFX = kurt_err_tlscb_PFX;
st_prtb_results.kurt_err_bezier_PFX = kurt_err_bezier_PFX;
st_prtb_results.kurt_err_lsqexp_APX = kurt_err_lsqexp_APX;
st_prtb_results.kurt_err_olscb_APX = kurt_err_olscb_APX;
st_prtb_results.kurt_err_tlscb_APX = kurt_err_tlscb_APX;
st_prtb_results.kurt_err_bezier_APX = kurt_err_bezier_APX;
st_prtb_results.kurt_err_lsqexp_APXF = kurt_err_lsqexp_APXF;
st_prtb_results.kurt_err_olscb_APXF = kurt_err_olscb_APXF;
st_prtb_results.kurt_err_tlscb_APXF = kurt_err_tlscb_APXF;
st_prtb_results.kurt_err_bezier_APXF = kurt_err_bezier_APXF;

if ( b_do_extr )
    st_prtb_results.mean_err_olscb_PF_extr = mean_err_olscb_PF_extr;
    st_prtb_results.mean_err_tlscb_PF_extr = mean_err_tlscb_PF_extr;
    st_prtb_results.mean_err_olscb_PF_extr2 = mean_err_olscb_PF_extr2;
    st_prtb_results.mean_err_olscb_PFX_extr = mean_err_olscb_PFX_extr;
    st_prtb_results.mean_err_tlscb_PFX_extr = mean_err_tlscb_PFX_extr;
%     st_prtb_results.mean_err_olscb_PFX_extr2 = mean_err_olscb_PFX_extr2;
    st_prtb_results.mean_err_olscb_APX_extr = mean_err_olscb_APX_extr;
    st_prtb_results.mean_err_tlscb_APX_extr = mean_err_tlscb_APX_extr;
    st_prtb_results.mean_err_olscb_APXF_extr = mean_err_olscb_APXF_extr;
    st_prtb_results.mean_err_tlscb_APXF_extr = mean_err_tlscb_APXF_extr;
    
    st_prtb_results.std_err_olscb_PF_extr = std_err_olscb_PF_extr;
    st_prtb_results.std_err_tlscb_PF_extr = std_err_tlscb_PF_extr;
    st_prtb_results.std_err_olscb_PF_extr2 = std_err_olscb_PF_extr2;
    st_prtb_results.std_err_olscb_PFX_extr = std_err_olscb_PFX_extr;
    st_prtb_results.std_err_tlscb_PFX_extr = std_err_tlscb_PFX_extr;
%     st_prtb_results.std_err_olscb_PFX_extr2 = std_err_olscb_PFX_extr2;
    st_prtb_results.std_err_olscb_APX_extr = std_err_olscb_APX_extr;
    st_prtb_results.std_err_tlscb_APX_extr = std_err_tlscb_APX_extr;
    st_prtb_results.std_err_olscb_APXF_extr = std_err_olscb_APXF_extr;
    st_prtb_results.std_err_tlscb_APXF_extr = std_err_tlscb_APXF_extr;
    
    st_prtb_results.skew_err_olscb_PF_extr = skew_err_olscb_PF_extr;
    st_prtb_results.skew_err_tlscb_PF_extr = skew_err_tlscb_PF_extr;
    st_prtb_results.skew_err_olscb_PF_extr2 = skew_err_olscb_PF_extr2;
    st_prtb_results.skew_err_olscb_PFX_extr = skew_err_olscb_PFX_extr;
    st_prtb_results.skew_err_tlscb_PFX_extr = skew_err_tlscb_PFX_extr;
%     st_prtb_results.skew_err_olscb_PFX_extr2 = skew_err_olscb_PFX_extr2;
    st_prtb_results.skew_err_olscb_APX_extr = skew_err_olscb_APX_extr;
    st_prtb_results.skew_err_tlscb_APX_extr = skew_err_tlscb_APX_extr;
    st_prtb_results.skew_err_olscb_APXF_extr = skew_err_olscb_APXF_extr;
    st_prtb_results.skew_err_tlscb_APXF_extr = skew_err_tlscb_APXF_extr;
    
    st_prtb_results.kurt_err_olscb_PF_extr = kurt_err_olscb_PF_extr;
    st_prtb_results.kurt_err_tlscb_PF_extr = kurt_err_tlscb_PF_extr;
    st_prtb_results.kurt_err_olscb_PF_extr2 = kurt_err_olscb_PF_extr2;
    st_prtb_results.kurt_err_olscb_PFX_extr = kurt_err_olscb_PFX_extr;
    st_prtb_results.kurt_err_tlscb_PFX_extr = kurt_err_tlscb_PFX_extr;
%     st_prtb_results.kurt_err_olscb_PFX_extr2 = kurt_err_olscb_PFX_extr2;
    st_prtb_results.kurt_err_olscb_APX_extr = kurt_err_olscb_APX_extr;
    st_prtb_results.kurt_err_tlscb_APX_extr = kurt_err_tlscb_APX_extr;
    st_prtb_results.kurt_err_olscb_APXF_extr = kurt_err_olscb_APXF_extr;
    st_prtb_results.kurt_err_tlscb_APXF_extr = kurt_err_tlscb_APXF_extr;
end


st_prtb_results.mean_nw_err_lsqexp_PF = mean_nw_err_lsqexp_PF;
st_prtb_results.mean_nw_err_olscb_PF = mean_nw_err_olscb_PF;
st_prtb_results.mean_nw_err_tlscb_PF = mean_nw_err_tlscb_PF;
st_prtb_results.mean_nw_err_bezier_PF = mean_nw_err_bezier_PF;
st_prtb_results.mean_nw_err_lsqexp_PFX = mean_nw_err_lsqexp_PFX;
st_prtb_results.mean_nw_err_olscb_PFX = mean_nw_err_olscb_PFX;
st_prtb_results.mean_nw_err_tlscb_PFX = mean_nw_err_tlscb_PFX;
st_prtb_results.mean_nw_err_bezier_PFX = mean_nw_err_bezier_PFX;
st_prtb_results.mean_nw_err_lsqexp_APX = mean_nw_err_lsqexp_APX;
st_prtb_results.mean_nw_err_olscb_APX = mean_nw_err_olscb_APX;
st_prtb_results.mean_nw_err_tlscb_APX = mean_nw_err_tlscb_APX;
st_prtb_results.mean_nw_err_bezier_APX = mean_nw_err_bezier_APX;
st_prtb_results.mean_nw_err_lsqexp_APXF = mean_nw_err_lsqexp_APXF;
st_prtb_results.mean_nw_err_olscb_APXF = mean_nw_err_olscb_APXF;
st_prtb_results.mean_nw_err_tlscb_APXF = mean_nw_err_tlscb_APXF;
st_prtb_results.mean_nw_err_bezier_APXF = mean_nw_err_bezier_APXF;

st_prtb_results.std_nw_err_lsqexp_PF = std_nw_err_lsqexp_PF;
st_prtb_results.std_nw_err_olscb_PF = std_nw_err_olscb_PF;
st_prtb_results.std_nw_err_tlscb_PF = std_nw_err_tlscb_PF;
st_prtb_results.std_nw_err_bezier_PF = std_nw_err_bezier_PF;
st_prtb_results.std_nw_err_lsqexp_PFX = std_nw_err_lsqexp_PFX;
st_prtb_results.std_nw_err_olscb_PFX = std_nw_err_olscb_PFX;
st_prtb_results.std_nw_err_tlscb_PFX = std_nw_err_tlscb_PFX;
st_prtb_results.std_nw_err_bezier_PFX = std_nw_err_bezier_PFX;
st_prtb_results.std_nw_err_lsqexp_APX = std_nw_err_lsqexp_APX;
st_prtb_results.std_nw_err_olscb_APX = std_nw_err_olscb_APX;
st_prtb_results.std_nw_err_tlscb_APX = std_nw_err_tlscb_APX;
st_prtb_results.std_nw_err_bezier_APX = std_nw_err_bezier_APX;
st_prtb_results.std_nw_err_lsqexp_APXF = std_nw_err_lsqexp_APXF;
st_prtb_results.std_nw_err_olscb_APXF = std_nw_err_olscb_APXF;
st_prtb_results.std_nw_err_tlscb_APXF = std_nw_err_tlscb_APXF;
st_prtb_results.std_nw_err_bezier_APXF = std_nw_err_bezier_APXF;

st_prtb_results.skew_nw_err_lsqexp_PF = skew_nw_err_lsqexp_PF;
st_prtb_results.skew_nw_err_olscb_PF = skew_nw_err_olscb_PF;
st_prtb_results.skew_nw_err_tlscb_PF = skew_nw_err_tlscb_PF;
st_prtb_results.skew_nw_err_bezier_PF = skew_nw_err_bezier_PF;
st_prtb_results.skew_nw_err_lsqexp_PFX = skew_nw_err_lsqexp_PFX;
st_prtb_results.skew_nw_err_olscb_PFX = skew_nw_err_olscb_PFX;
st_prtb_results.skew_nw_err_tlscb_PFX = skew_nw_err_tlscb_PFX;
st_prtb_results.skew_nw_err_bezier_PFX = skew_nw_err_bezier_PFX;
st_prtb_results.skew_nw_err_lsqexp_APX = skew_nw_err_lsqexp_APX;
st_prtb_results.skew_nw_err_olscb_APX = skew_nw_err_olscb_APX;
st_prtb_results.skew_nw_err_tlscb_APX = skew_nw_err_tlscb_APX;
st_prtb_results.skew_nw_err_bezier_APX = skew_nw_err_bezier_APX;
st_prtb_results.skew_nw_err_lsqexp_APXF = skew_nw_err_lsqexp_APXF;
st_prtb_results.skew_nw_err_olscb_APXF = skew_nw_err_olscb_APXF;
st_prtb_results.skew_nw_err_tlscb_APXF = skew_nw_err_tlscb_APXF;
st_prtb_results.skew_nw_err_bezier_APXF = skew_nw_err_bezier_APXF;

st_prtb_results.kurt_nw_err_lsqexp_PF = kurt_nw_err_lsqexp_PF;
st_prtb_results.kurt_nw_err_olscb_PF = kurt_nw_err_olscb_PF;
st_prtb_results.kurt_nw_err_tlscb_PF = kurt_nw_err_tlscb_PF;
st_prtb_results.kurt_nw_err_bezier_PF = kurt_nw_err_bezier_PF;
st_prtb_results.kurt_nw_err_lsqexp_PFX = kurt_nw_err_lsqexp_PFX;
st_prtb_results.kurt_nw_err_olscb_PFX = kurt_nw_err_olscb_PFX;
st_prtb_results.kurt_nw_err_tlscb_PFX = kurt_nw_err_tlscb_PFX;
st_prtb_results.kurt_nw_err_bezier_PFX = kurt_nw_err_bezier_PFX;
st_prtb_results.kurt_nw_err_lsqexp_APX = kurt_nw_err_lsqexp_APX;
st_prtb_results.kurt_nw_err_olscb_APX = kurt_nw_err_olscb_APX;
st_prtb_results.kurt_nw_err_tlscb_APX = kurt_nw_err_tlscb_APX;
st_prtb_results.kurt_nw_err_bezier_APX = kurt_nw_err_bezier_APX;
st_prtb_results.kurt_nw_err_lsqexp_APXF = kurt_nw_err_lsqexp_APXF;
st_prtb_results.kurt_nw_err_olscb_APXF = kurt_nw_err_olscb_APXF;
st_prtb_results.kurt_nw_err_tlscb_APXF = kurt_nw_err_tlscb_APXF;
st_prtb_results.kurt_nw_err_bezier_APXF = kurt_nw_err_bezier_APXF;

if ( b_do_extr )
    st_prtb_results.mean_nw_err_olscb_PF_extr = mean_nw_err_olscb_PF_extr;
    st_prtb_results.mean_nw_err_tlscb_PF_extr = mean_nw_err_tlscb_PF_extr;
    st_prtb_results.mean_nw_err_olscb_PF_extr2 = mean_nw_err_olscb_PF_extr2;
    st_prtb_results.mean_nw_err_olscb_PFX_extr = mean_nw_err_olscb_PFX_extr;
    st_prtb_results.mean_nw_err_tlscb_PFX_extr = mean_nw_err_tlscb_PFX_extr;
%     st_prtb_results.mean_nw_err_olscb_PFX_extr2 = mean_nw_err_olscb_PFX_extr2;   
    st_prtb_results.mean_nw_err_olscb_APX_extr = mean_nw_err_olscb_APX_extr;
    st_prtb_results.mean_nw_err_tlscb_APX_extr = mean_nw_err_tlscb_APX_extr;
    st_prtb_results.mean_nw_err_olscb_APXF_extr = mean_nw_err_olscb_APXF_extr;
    st_prtb_results.mean_nw_err_tlscb_APXF_extr = mean_nw_err_tlscb_APXF_extr;
    
    st_prtb_results.std_nw_err_olscb_PF_extr = std_nw_err_olscb_PF_extr;
    st_prtb_results.std_nw_err_tlscb_PF_extr = std_nw_err_tlscb_PF_extr;
    st_prtb_results.std_nw_err_olscb_PF_extr2 = std_nw_err_olscb_PF_extr2;
    st_prtb_results.std_nw_err_olscb_PFX_extr = std_nw_err_olscb_PFX_extr;
    st_prtb_results.std_nw_err_tlscb_PFX_extr = std_nw_err_tlscb_PFX_extr;
%     st_prtb_results.std_nw_err_olscb_PFX_extr2 = std_nw_err_olscb_PFX_extr2;
    st_prtb_results.std_nw_err_olscb_APX_extr = std_nw_err_olscb_APX_extr;
    st_prtb_results.std_nw_err_tlscb_APX_extr = std_nw_err_tlscb_APX_extr;
    st_prtb_results.std_nw_err_olscb_APXF_extr = std_nw_err_olscb_APXF_extr;
    st_prtb_results.std_nw_err_tlscb_APXF_extr = std_nw_err_tlscb_APXF_extr;
    
    st_prtb_results.skew_nw_err_olscb_PF_extr = skew_nw_err_olscb_PF_extr;
    st_prtb_results.skew_nw_err_tlscb_PF_extr = skew_nw_err_tlscb_PF_extr;
    st_prtb_results.skew_nw_err_olscb_PF_extr2 = skew_nw_err_olscb_PF_extr2;
    st_prtb_results.skew_nw_err_olscb_PFX_extr = skew_nw_err_olscb_PFX_extr;
    st_prtb_results.skew_nw_err_tlscb_PFX_extr = skew_nw_err_tlscb_PFX_extr;
%     st_prtb_results.skew_nw_err_olscb_PFX_extr2 = skew_nw_err_olscb_PFX_extr2;
    st_prtb_results.skew_nw_err_olscb_APX_extr = skew_nw_err_olscb_APX_extr;
    st_prtb_results.skew_nw_err_tlscb_APX_extr = skew_nw_err_tlscb_APX_extr;
    st_prtb_results.skew_nw_err_olscb_APXF_extr = skew_nw_err_olscb_APXF_extr;
    st_prtb_results.skew_nw_err_tlscb_APXF_extr = skew_nw_err_tlscb_APXF_extr;
    
    st_prtb_results.kurt_nw_err_olscb_PF_extr = kurt_nw_err_olscb_PF_extr;
    st_prtb_results.kurt_nw_err_tlscb_PF_extr = kurt_nw_err_tlscb_PF_extr;
    st_prtb_results.kurt_nw_err_olscb_PF_extr2 = kurt_nw_err_olscb_PF_extr2;
    st_prtb_results.kurt_nw_err_olscb_PFX_extr = kurt_nw_err_olscb_PFX_extr;
    st_prtb_results.kurt_nw_err_tlscb_PFX_extr = kurt_nw_err_tlscb_PFX_extr;
%     st_prtb_results.kurt_nw_err_olscb_PFX_extr2 = kurt_nw_err_olscb_PFX_extr2;
    st_prtb_results.kurt_nw_err_olscb_APX_extr = kurt_nw_err_olscb_APX_extr;
    st_prtb_results.kurt_nw_err_tlscb_APX_extr = kurt_nw_err_tlscb_APX_extr;
    st_prtb_results.kurt_nw_err_olscb_APXF_extr = kurt_nw_err_olscb_APXF_extr;
    st_prtb_results.kurt_nw_err_tlscb_APXF_extr = kurt_nw_err_tlscb_APXF_extr;
end


st_prtb_results.v_c_nw_err_lsqexp_PF = v_c_nw_err_lsqexp_PF;
st_prtb_results.v_c_nw_err_olscb_PF = v_c_nw_err_olscb_PF;
st_prtb_results.v_c_nw_err_tlscb_PF = v_c_nw_err_tlscb_PF;
st_prtb_results.v_c_nw_err_bezier_PF = v_c_nw_err_bezier_PF;
st_prtb_results.v_c_nw_err_lsqexp_PFX = v_c_nw_err_lsqexp_PFX;
st_prtb_results.v_c_nw_err_olscb_PFX = v_c_nw_err_olscb_PFX;
st_prtb_results.v_c_nw_err_tlscb_PFX = v_c_nw_err_tlscb_PFX;
st_prtb_results.v_c_nw_err_bezier_PFX = v_c_nw_err_bezier_PFX;
st_prtb_results.v_c_nw_err_lsqexp_APX = v_c_nw_err_lsqexp_APX;
st_prtb_results.v_c_nw_err_olscb_APX = v_c_nw_err_olscb_APX;
st_prtb_results.v_c_nw_err_tlscb_APX = v_c_nw_err_tlscb_APX;
st_prtb_results.v_c_nw_err_bezier_APX = v_c_nw_err_bezier_APX;
st_prtb_results.v_c_nw_err_lsqexp_APXF = v_c_nw_err_lsqexp_APXF;
st_prtb_results.v_c_nw_err_olscb_APXF = v_c_nw_err_olscb_APXF;
st_prtb_results.v_c_nw_err_tlscb_APXF = v_c_nw_err_tlscb_APXF;
st_prtb_results.v_c_nw_err_bezier_APXF = v_c_nw_err_bezier_APXF;
if ( b_do_extr )
    st_prtb_results.v_c_nw_err_olscb_PF_extr = v_c_nw_err_olscb_PF_extr;
    st_prtb_results.v_c_nw_err_tlscb_PF_extr = v_c_nw_err_tlscb_PF_extr;
    st_prtb_results.v_c_nw_err_olscb_PF_extr2 = v_c_nw_err_olscb_PF_extr2;
    st_prtb_results.v_c_nw_err_olscb_PFX_extr = v_c_nw_err_olscb_PFX_extr;
    st_prtb_results.v_c_nw_err_tlscb_PFX_extr = v_c_nw_err_tlscb_PFX_extr;
%     st_prtb_results.v_c_nw_err_olscb_PFX_extr2 = v_c_nw_err_olscb_PFX_extr2;
    st_prtb_results.v_c_nw_err_olscb_APX_extr = v_c_nw_err_olscb_APX_extr;
    st_prtb_results.v_c_nw_err_tlscb_APX_extr = v_c_nw_err_tlscb_APX_extr;
    st_prtb_results.v_c_nw_err_olscb_APXF_extr = v_c_nw_err_olscb_APXF_extr;
    st_prtb_results.v_c_nw_err_tlscb_APXF_extr = v_c_nw_err_tlscb_APXF_extr;
end


st_prtb_results.a_prctl_mean_err_lsqexp_PF = a_prctl_mean_err_lsqexp_PF;
st_prtb_results.a_prctl_mean_err_olscb_PF = a_prctl_mean_err_olscb_PF;
st_prtb_results.a_prctl_mean_err_tlscb_PF = a_prctl_mean_err_tlscb_PF;
st_prtb_results.a_prctl_mean_err_bezier_PF = a_prctl_mean_err_bezier_PF;
st_prtb_results.a_prctl_mean_err_lsqexp_PFX = a_prctl_mean_err_lsqexp_PFX;
st_prtb_results.a_prctl_mean_err_olscb_PFX = a_prctl_mean_err_olscb_PFX;
st_prtb_results.a_prctl_mean_err_tlscb_PFX = a_prctl_mean_err_tlscb_PFX;
st_prtb_results.a_prctl_mean_err_bezier_PFX = a_prctl_mean_err_bezier_PFX;
st_prtb_results.a_prctl_mean_err_lsqexp_APX = a_prctl_mean_err_lsqexp_APX;
st_prtb_results.a_prctl_mean_err_olscb_APX = a_prctl_mean_err_olscb_APX;
st_prtb_results.a_prctl_mean_err_tlscb_APX = a_prctl_mean_err_tlscb_APX;
st_prtb_results.a_prctl_mean_err_bezier_APX = a_prctl_mean_err_bezier_APX;
st_prtb_results.a_prctl_mean_err_lsqexp_APXF = a_prctl_mean_err_lsqexp_APXF;
st_prtb_results.a_prctl_mean_err_olscb_APXF = a_prctl_mean_err_olscb_APXF;
st_prtb_results.a_prctl_mean_err_tlscb_APXF = a_prctl_mean_err_tlscb_APXF;
st_prtb_results.a_prctl_mean_err_bezier_APXF = a_prctl_mean_err_bezier_APXF;
if ( b_do_extr )    
    st_prtb_results.a_prctl_mean_err_olscb_PF_extr = a_prctl_mean_err_olscb_PF_extr;
    st_prtb_results.a_prctl_mean_err_tlscb_PF_extr = a_prctl_mean_err_tlscb_PF_extr;
    st_prtb_results.a_prctl_mean_err_olscb_PF_extr2 = a_prctl_mean_err_olscb_PF_extr2;    
    st_prtb_results.a_prctl_mean_err_olscb_PFX_extr = a_prctl_mean_err_olscb_PFX_extr;
    st_prtb_results.a_prctl_mean_err_tlscb_PFX_extr = a_prctl_mean_err_tlscb_PFX_extr;
%     st_prtb_results.a_prctl_mean_err_olscb_PFX_extr2 = a_prctl_mean_err_olscb_PFX_extr2;
    st_prtb_results.a_prctl_mean_err_olscb_APX_extr = a_prctl_mean_err_olscb_APX_extr;
    st_prtb_results.a_prctl_mean_err_tlscb_APX_extr = a_prctl_mean_err_tlscb_APX_extr;
    st_prtb_results.a_prctl_mean_err_olscb_APXF_extr = a_prctl_mean_err_olscb_APXF_extr;
    st_prtb_results.a_prctl_mean_err_tlscb_APXF_extr = a_prctl_mean_err_tlscb_APXF_extr;
end


st_prtb_results.v_c_lsqexp_PF_min = v_c_lsqexp_PF_min;
st_prtb_results.v_c_lsqexp_PF_max = v_c_lsqexp_PF_max;
st_prtb_results.v_c_olscb_PF_min = v_c_olscb_PF_min;
st_prtb_results.v_c_olscb_PF_max = v_c_olscb_PF_max;
st_prtb_results.v_c_tlscb_PF_min = v_c_tlscb_PF_min;
st_prtb_results.v_c_tlscb_PF_max = v_c_tlscb_PF_max;
st_prtb_results.v_c_bezier_PF_min = v_c_bezier_PF_min;
st_prtb_results.v_c_bezier_PF_max = v_c_bezier_PF_max;
st_prtb_results.v_c_lsqexp_PFX_min = v_c_lsqexp_PFX_min;
st_prtb_results.v_c_lsqexp_PFX_max = v_c_lsqexp_PFX_max;
st_prtb_results.v_c_olscb_PFX_min = v_c_olscb_PFX_min;
st_prtb_results.v_c_olscb_PFX_max = v_c_olscb_PFX_max;
st_prtb_results.v_c_tlscb_PFX_min = v_c_tlscb_PFX_min;
st_prtb_results.v_c_tlscb_PFX_max = v_c_tlscb_PFX_max;
st_prtb_results.v_c_bezier_PFX_min = v_c_bezier_PFX_min;
st_prtb_results.v_c_bezier_PFX_max = v_c_bezier_PFX_max;
st_prtb_results.v_c_lsqexp_APX_min = v_c_lsqexp_APX_min;
st_prtb_results.v_c_lsqexp_APX_max = v_c_lsqexp_APX_max;
st_prtb_results.v_c_olscb_APX_min = v_c_olscb_APX_min;
st_prtb_results.v_c_olscb_APX_max = v_c_olscb_APX_max;
st_prtb_results.v_c_tlscb_APX_min = v_c_tlscb_APX_min;
st_prtb_results.v_c_tlscb_APX_max = v_c_tlscb_APX_max;
st_prtb_results.v_c_bezier_APX_min = v_c_bezier_APX_min;
st_prtb_results.v_c_bezier_APX_max = v_c_bezier_APX_max;
st_prtb_results.v_c_lsqexp_APXF_min = v_c_lsqexp_APXF_min;
st_prtb_results.v_c_lsqexp_APXF_max = v_c_lsqexp_APXF_max;
st_prtb_results.v_c_olscb_APXF_min = v_c_olscb_APXF_min;
st_prtb_results.v_c_olscb_APXF_max = v_c_olscb_APXF_max;
st_prtb_results.v_c_tlscb_APXF_min = v_c_tlscb_APXF_min;
st_prtb_results.v_c_tlscb_APXF_max = v_c_tlscb_APXF_max;
st_prtb_results.v_c_bezier_APXF_min = v_c_bezier_APXF_min;
st_prtb_results.v_c_bezier_APXF_max = v_c_bezier_APXF_max;
if ( b_do_extr )
    st_prtb_results.v_c_olscb_PF_extr_min = v_c_olscb_PF_extr_min;
    st_prtb_results.v_c_olscb_PF_extr_max = v_c_olscb_PF_extr_max;
    st_prtb_results.v_c_tlscb_PF_extr_min = v_c_tlscb_PF_extr_min;
    st_prtb_results.v_c_tlscb_PF_extr_max = v_c_tlscb_PF_extr_max;
    st_prtb_results.v_c_olscb_PF_extr2_min = v_c_olscb_PF_extr2_min;
    st_prtb_results.v_c_olscb_PF_extr2_max = v_c_olscb_PF_extr2_max;
    st_prtb_results.v_c_olscb_PFX_extr_min = v_c_olscb_PFX_extr_min;
    st_prtb_results.v_c_olscb_PFX_extr_max = v_c_olscb_PFX_extr_max;
    st_prtb_results.v_c_tlscb_PFX_extr_min = v_c_tlscb_PFX_extr_min;
    st_prtb_results.v_c_tlscb_PFX_extr_max = v_c_tlscb_PFX_extr_max;
%     st_prtb_results.v_c_olscb_PFX_extr2_min = v_c_olscb_PFX_extr2_min;
%     st_prtb_results.v_c_olscb_PFX_extr2_max = v_c_olscb_PFX_extr2_max;
    st_prtb_results.v_c_olscb_APX_extr_min = v_c_olscb_APX_extr_min;
    st_prtb_results.v_c_olscb_APX_extr_max = v_c_olscb_APX_extr_max;
    st_prtb_results.v_c_tlscb_APX_extr_min = v_c_tlscb_APX_extr_min;
    st_prtb_results.v_c_tlscb_APX_extr_max = v_c_tlscb_APX_extr_max;
    st_prtb_results.v_c_olscb_APXF_extr_min = v_c_olscb_APXF_extr_min;
    st_prtb_results.v_c_olscb_APXF_extr_max = v_c_olscb_APXF_extr_max;
    st_prtb_results.v_c_tlscb_APXF_extr_min = v_c_tlscb_APXF_extr_min;
    st_prtb_results.v_c_tlscb_APXF_extr_max = v_c_tlscb_APXF_extr_max;
end


st_prtb_results.err_min_lsqexp_PF = ( v_c_lsqexp_PF_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_lsqexp_PF = ( v_c_lsqexp_PF_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_olscb_PF = ( v_c_olscb_PF_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_olscb_PF = ( v_c_olscb_PF_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_tlscb_PF = ( v_c_tlscb_PF_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_tlscb_PF = ( v_c_tlscb_PF_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_bezier_PF = ( v_c_bezier_PF_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_bezier_PF = ( v_c_bezier_PF_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_lsqexp_PFX = ( v_c_lsqexp_PFX_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_lsqexp_PFX = ( v_c_lsqexp_PFX_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_olscb_PFX = ( v_c_olscb_PFX_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_olscb_PFX = ( v_c_olscb_PFX_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_tlscb_PFX = ( v_c_tlscb_PFX_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_tlscb_PFX = ( v_c_tlscb_PFX_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_bezier_PFX = ( v_c_bezier_PFX_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_bezier_PFX = ( v_c_bezier_PFX_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_lsqexp_APX = ( v_c_lsqexp_APX_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_lsqexp_APX = ( v_c_lsqexp_APX_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_olscb_APX = ( v_c_olscb_APX_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_olscb_APX = ( v_c_olscb_APX_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_tlscb_APX = ( v_c_tlscb_APX_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_tlscb_APX = ( v_c_tlscb_APX_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_bezier_APX = ( v_c_bezier_APX_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_bezier_APX = ( v_c_bezier_APX_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_lsqexp_APXF = ( v_c_lsqexp_APXF_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_lsqexp_APXF = ( v_c_lsqexp_APXF_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_olscb_APXF = ( v_c_olscb_APXF_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_olscb_APXF = ( v_c_olscb_APXF_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_tlscb_APXF = ( v_c_tlscb_APXF_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_tlscb_APXF = ( v_c_tlscb_APXF_max - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_min_bezier_APXF = ( v_c_bezier_APXF_min - v_c_ref ) ./ v_c_ref;
st_prtb_results.err_max_bezier_APXF = ( v_c_bezier_APXF_max - v_c_ref ) ./ v_c_ref;
if ( b_do_extr )
    st_prtb_results.err_min_olscb_PF_extr = ( v_c_olscb_PF_extr_min - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_max_olscb_PF_extr = ( v_c_olscb_PF_extr_max - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_min_tlscb_PF_extr = ( v_c_tlscb_PF_extr_min - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_max_tlscb_PF_extr = ( v_c_tlscb_PF_extr_max - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_min_olscb_PF_extr2 = ( v_c_olscb_PF_extr2_min - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_max_olscb_PF_extr2 = ( v_c_olscb_PF_extr2_max - v_c_ref ) ./ v_c_ref;    
    st_prtb_results.err_min_olscb_PFX_extr = ( v_c_olscb_PFX_extr_min - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_max_olscb_PFX_extr = ( v_c_olscb_PFX_extr_max - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_min_tlscb_PFX_extr = ( v_c_tlscb_PFX_extr_min - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_max_tlscb_PFX_extr = ( v_c_tlscb_PFX_extr_max - v_c_ref ) ./ v_c_ref;
%     st_prtb_results.err_min_olscb_PFX_extr2 = ( v_c_olscb_PFX_extr2_min - v_c_ref ) ./ v_c_ref;
%     st_prtb_results.err_max_olscb_PFX_extr2 = ( v_c_olscb_PFX_extr2_max - v_c_ref ) ./ v_c_ref;    
    st_prtb_results.err_min_olscb_APX_extr = ( v_c_olscb_APX_extr_min - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_max_olscb_APX_extr = ( v_c_olscb_APX_extr_max - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_min_tlscb_APX_extr = ( v_c_tlscb_APX_extr_min - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_max_tlscb_APX_extr = ( v_c_tlscb_APX_extr_max - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_min_olscb_APXF_extr = ( v_c_olscb_APXF_extr_min - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_max_olscb_APXF_extr = ( v_c_olscb_APXF_extr_max - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_min_tlscb_APXF_extr = ( v_c_tlscb_APXF_extr_min - v_c_ref ) ./ v_c_ref;
    st_prtb_results.err_max_tlscb_APXF_extr = ( v_c_tlscb_APXF_extr_max - v_c_ref ) ./ v_c_ref;
end


% Save stat results
st_prtb_results.v_k = v_k;
st_prtb_results.st_p = st_p;
st_prtb_results.Nitr = Nitr;
st_prtb_results.st_fn_shear_exact = st_fn_shear_exact;
st_prtb_results.v_zs = v_zs;
st_prtb_results.v_zs_err_width = v_zs_err_width;
st_prtb_results.Nz = Nz;
st_prtb_results.st_Dn = st_Dn;
st_prtb_results.v_c_ref = v_c_ref;
st_prtb_results.v_percentile = v_percentile;
st_prtb_results.itr_done_pf = itr_done_pf;
st_prtb_results.itr_done_apx = itr_done_apx;
save( s_filename, 'st_prtb_results' );


% Save backup
save( s_backup_filename, 'ca_backup' );


% figure(1);
% plot( v_k, mean_err_olscb, 'b', v_k, std_err_olscb, 'g--', v_k, skew_err_olscb, 'm-.', v_k, kurt_err_olscb, 'r:' );
% 
% figure(2);
% plot( v_k, mean_err_tlscb, 'b', v_k, std_err_tlscb, 'g--', v_k, skew_err_tlscb, 'm-.', v_k, kurt_err_tlscb, 'r:' );
% 
% figure(3);
% plot( v_k, mean_err_olscb_extr, 'b', v_k, std_err_olscb_extr, 'g--', v_k, skew_err_olscb_extr, 'm-.', v_k, kurt_err_olscb_extr, 'r:' );
% 
% figure(4);
% plot( v_k, mean_err_tlscb_extr, 'b', v_k, std_err_tlscb_extr, 'g--', v_k, skew_err_tlscb_extr, 'm-.', v_k, kurt_err_tlscb_extr, 'r:' );
% 
% 
% 
% norm( mean_err_olscb, 2 ) / numel( v_k )
% norm( mean_err_tlscb, 2 ) / numel( v_k )
% norm( mean_err_olscb_extr, 2 ) / numel( v_k )
% norm( mean_err_tlscb_extr, 2 ) / numel( v_k )

% 
%     % Find min and max lsqexp
%     v_lsqexp_min_idxf = v_c_lsqexp_min > st_shear_prtb_single.v_c_pf_lsqexp;
%     v_c_lsqexp_min( v_lsqexp_min_idxf ) =  st_shear_prtb_single.v_c_pf_lsqexp( v_lsqexp_min_idxf );
%     v_lsqexp_max_idxf = v_c_lsqexp_max < st_shear_prtb_single.v_c_pf_lsqexp;
%     v_c_lsqexp_max( v_lsqexp_max_idxf ) =  st_shear_prtb_single.v_c_pf_lsqexp( v_lsqexp_max_idxf );
% 
%     v_olscb_min_idxf = v_c_olscb_min > st_shear_prtb_single.v_c_pf_olscb;
%     v_c_olscb_min( v_olscb_min_idxf ) =  st_shear_prtb_single.v_c_pf_olscb( v_olscb_min_idxf );
%     v_olscb_max_idxf = v_c_olscb_max < st_shear_prtb_single.v_c_pf_olscb;
%     v_c_olscb_max( v_olscb_max_idxf ) =  st_shear_prtb_single.v_c_pf_olscb( v_olscb_max_idxf );
%     
%     v_tlscb_min_idxf = v_c_tlscb_min > st_shear_prtb_single.v_c_pf_tlscb;
%     v_c_tlscb_min( v_tlscb_min_idxf ) =  st_shear_prtb_single.v_c_pf_tlscb( v_tlscb_min_idxf );
%     v_tlscb_max_idxf = v_c_tlscb_max < st_shear_prtb_single.v_c_pf_tlscb;
%     v_c_tlscb_max( v_tlscb_max_idxf ) =  st_shear_prtb_single.v_c_pf_tlscb( v_tlscb_max_idxf );
% 
%     if ( b_do_extr )
%         v_olscb_extr_min_idxf = v_c_olscb_extr_min > st_shear_prtb_single.v_c_pf_olscb_extr;
%         v_c_olscb_extr_min( v_olscb_extr_min_idxf ) =  st_shear_prtb_single.v_c_pf_olscb_extr( v_olscb_extr_min_idxf );
%         v_olscb_extr_max_idxf = v_c_olscb_extr_max < st_shear_prtb_single.v_c_pf_olscb_extr;
%         v_c_olscb_extr_max( v_olscb_extr_max_idxf ) =  st_shear_prtb_single.v_c_pf_olscb_extr( v_olscb_extr_max_idxf );
% 
%         v_tlscb_extr_min_idxf = v_c_tlscb_extr_min > st_shear_prtb_single.v_c_pf_tlscb_extr;
%         v_c_tlscb_extr_min( v_tlscb_extr_min_idxf ) =  st_shear_prtb_single.v_c_pf_tlscb_extr( v_tlscb_extr_min_idxf );
%         v_tlscb_extr_max_idxf = v_c_tlscb_extr_max < st_shear_prtb_single.v_c_pf_tlscb_extr;
%         v_c_tlscb_extr_max( v_tlscb_extr_max_idxf ) =  st_shear_prtb_single.v_c_pf_tlscb_extr( v_tlscb_extr_max_idxf );
%         
%         v_olscb_extr2_min_idxf = v_c_olscb_extr2_min > st_shear_prtb_single.v_c_pf_olscb_extr2;
%         v_c_olscb_extr2_min( v_olscb_extr2_min_idxf ) =  st_shear_prtb_single.v_c_pf_olscb_extr2( v_olscb_extr2_min_idxf );
%         v_olscb_extr2_max_idxf = v_c_olscb_extr2_max < st_shear_prtb_single.v_c_pf_olscb_extr2;
%         v_c_olscb_extr2_max( v_olscb_extr2_max_idxf ) =  st_shear_prtb_single.v_c_pf_olscb_extr2( v_olscb_extr2_max_idxf );
%     end


end
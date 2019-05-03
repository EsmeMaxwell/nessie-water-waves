function [ st_prtb_single, st_err ] = fn_ww__sim_shear_prtb__do_single( st_Dn, st_res, v_zs, v_zs_exact_U, v_zs_err_width, v_k, b_do_extr, st_p )
%fn_ww__sim_shear_prtb__do_single: Sim (parfor fn) single run for err-in-shear problem
%
%   [ st_shear_prtb_single, st_err ] = fn_ww__sim_shear_prtb__do_single( st_Dn, st_res, v_zs, v_zs_exact_U, v_zs_err_width, v_k, b_do_extr, st_p )
%
% Performs a single run for the err-in-shear-profile simulation. ONLY call
% from fn_ww__sim_shear_prtb__do_runs().
% 
% See also
%   fn_ww__sim_shear_prtb__do_runs()



set(0,'defaulttextinterpreter','latex');

b_show_warnings = true;

upperlim_relerr = 3.0;
upperlim_clpf_apx_diff = 20;
Nextrpts = 4;
hs = -v_zs(1);

if ( nargout > 1 )
    st_err = struct;
    st_err.b_ok = true;
end


%--------------------------------------------------------------------------
% Create perturbed discrete shear profile
%--------------------------------------------------------------------------

% First we create the perturbed shear "measurement"
v_zs_candidate_U = v_zs_exact_U + v_zs_err_width .* randn( size( v_zs_exact_U ) );


%--------------------------------------------------------------------------
% Do the fits
%--------------------------------------------------------------------------

% Do LSQ MATLAB exp fit
v_lsqexp_fit_param = lsqcurvefit( st_res.fn__exp_fit, [1 1 0], v_zs, v_zs_candidate_U, [ ], [ ],  optimoptions( 'lsqcurvefit', 'MaxFunctionEvaluations', 1500, 'Display', 'off' ) );

% Do the OLS fits
[ st_fit_olscb_PF ] = fn_ww__calc_fit__ols( st_res.st_fit_cb_interp_PF, st_res.st_fit_olscb_mtrx, v_zs_candidate_U, st_p.h - hs );
[ st_fit_olscb_APX ] = fn_ww__calc_fit__ols( st_res.st_fit_cb_interp_APX, st_res.st_fit_olscb_mtrx, v_zs_candidate_U, st_p.h - hs );

% Do the TLS fit
[ st_fit_tlscb_PF ] = fn_ww__calc_fit__tls( st_res.st_fit_cb_interp_PF, st_res.st_fit_tlscb_mtrx, v_zs_candidate_U, st_p.h - hs );
[ st_fit_tlscb_APX ] = fn_ww__calc_fit__tls( st_res.st_fit_cb_interp_APX, st_res.st_fit_tlscb_mtrx, v_zs_candidate_U, st_p.h - hs );

% Do the Bezier fit (no extrapolation... the b_extr test below is somewhat
% redunant for the Bezier stuff but no matter)
if ( ~b_do_extr )
    [ st_fit_bezier_PF ] = fn_ww__calc_fit__bezier( st_Dn.v_zm, v_zs, v_zs_candidate_U, st_p.h, 0, 0 );
    [ st_fit_bezier_APX ] = fn_ww__calc_fit__bezier( st_res.v_z_cc, v_zs, v_zs_candidate_U, st_p.h, 0, 0 );
end


if ( b_do_extr )
    
    % Do the extrapolated OLS fit
    [ v_zs_candidate_U_extr ] = fn_ww__calc_fit__prep_lin_surf_extrapolate( v_zs, v_zs_candidate_U, Nextrpts );
    [ st_fit_olscb_PF_extr ] = fn_ww__calc_fit__ols( st_res.st_fit_cb_interp_PF_extr, st_res.st_fit_olscb_mtrx_extr, v_zs_candidate_U_extr, st_p.h );
    [ st_fit_olscb_APX_extr ] = fn_ww__calc_fit__ols( st_res.st_fit_cb_interp_APX_extr, st_res.st_fit_olscb_mtrx_extr, v_zs_candidate_U_extr, st_p.h );
    
    % Do the extrapolated TLS fit
    [ v_zs_candidate_U_extr ] = fn_ww__calc_fit__prep_lin_surf_extrapolate( v_zs, v_zs_candidate_U, Nextrpts );
    [ st_fit_tlscb_PF_extr ] = fn_ww__calc_fit__tls( st_res.st_fit_cb_interp_PF_extr, st_res.st_fit_tlscb_mtrx_extr, v_zs_candidate_U_extr, st_p.h );
    [ st_fit_tlscb_APX_extr ] = fn_ww__calc_fit__tls( st_res.st_fit_cb_interp_APX_extr, st_res.st_fit_tlscb_mtrx_extr, v_zs_candidate_U_extr, st_p.h );

    % Do the Bezier fit
    [ st_fit_bezier_PF ] = fn_ww__calc_fit__bezier( st_Dn.v_zm, v_zs, v_zs_candidate_U, st_p.h, hs, Nextrpts );
    [ st_fit_bezier_APX ] = fn_ww__calc_fit__bezier( st_res.v_z_cc, v_zs, v_zs_candidate_U, st_p.h, hs, Nextrpts );
    
%     plot( st_Dn.v_zm, st_fit_bezier_PF.v_y, 'b' );
%     hold on;
%     plot( st_res.v_z_cc, st_fit_bezier_APX.v_y, 'm--' );
%     hold off;
%     error();
    
    
    % We don't tend to use extr2, so I'll not code the APX stuff for it,
    % only PF fit.
    
    % Do OLS extrapolate type II (do the OLS fit but interpolate using the
    % original CB basis of measurements points so we can determine the
    % derivative at h_s. We then readjust the points above h_s in the
    % original OLS fit to be a linear extrapolation.    
    [ st_fit_olscb_minus_h_to_hs_domain ] = fn_ww__calc_fit__ols( st_res.st_fit_cb, st_res.st_fit_olscb_mtrx, v_zs_candidate_U, st_p.h - hs );
    v_hs_layer_idxs = find( st_res.st_fit_cb_interp_PF.v_zs >= -hs );  % Note, this is on a different domain -- the Chebyshev domain from [-h, 0]; the fit domain was [-h,-h_s]!!
    v_hs_layer_idxs = [ v_hs_layer_idxs; v_hs_layer_idxs(end)+1 ];  % Need to drop one more to make sure we actually join at the correct place
    m_hs = st_fit_olscb_minus_h_to_hs_domain.v_zs_fit_dU(1);
    val_hs = st_fit_olscb_minus_h_to_hs_domain.v_zs_fit_U(1);
    loc_hs = st_res.st_fit_cb_interp_PF.v_zs( v_hs_layer_idxs(end) );  % Note this is almost hs but not quite, will be slightly lower  
    v_lin_extr = m_hs * ( st_res.st_fit_cb_interp_PF.v_zs( v_hs_layer_idxs ) - loc_hs ) + val_hs;
    
%     figure(1);
%     plot( st_Dn.v_zm, st_res.st_fn_shear_exact.fn_U( st_Dn.v_zm ), 'b' );    
%     hold on;
%     plot( v_zs, st_fit_olscb_minus_h_to_hs_domain.v_zs_fit_U, 'r:', 'LineWidth', 1.5 );
%     plot( st_res.st_fit_cb_interp_PF.v_zs, 0 .* st_res.st_fit_cb_interp_PF.v_zs, 'k*' );  % This shoudl be chebyshev on [-h,0]
%     plot( v_zs, 0.1 .* ones( size( v_zs ) ), 'm*' );  % This shoudl be chebyshev on [-h,-h_s]
%     plot( st_res.st_fit_cb_interp_PF.v_zs( v_hs_layer_idxs ), v_lin_extr, 'g--', 'LineWidth', 1.5 ); % The extrapolation, hopefully
%     hold off; 
%     
%     figure(2)
%     plot( st_Dn.v_zm, st_res.st_fn_shear_exact.fn_dU( st_Dn.v_zm ), 'b' );    
%     hold on;
%     plot( v_zs, st_fit_olscb_minus_h_to_hs_domain.v_zs_fit_dU, 'r:', 'LineWidth', 1.5 );
%     plot( st_res.st_fit_cb_interp_PF.v_zs, 0 .* st_res.st_fit_cb_interp_PF.v_zs, 'k*' );  % This shoudl be chebyshev on [-h,0]
%     plot( v_zs, 5 .* ones( size( v_zs ) ), 'm*' );  % This shoudl be chebyshev on [-h,-h_s]
%     plot( st_res.st_fit_cb_interp_PF.v_zs( v_hs_layer_idxs ), m_hs * ones( size( v_hs_layer_idxs ) ), 'g--', 'LineWidth', 1.5 ); % The extrapolation, hopefully
%     hold off;    
%     error()    
end

% NOTE!  Need to tidy-up/fix all the calls using st_Dn to allow simple
% vector!

% Define vector struct for Curtis-Clenshaw vector
st_vec_cc = fn_ww__setup__create_st_vec( st_res.v_z_cc );

% For clarity, and that's the only benefit, we take the vector components
% out of st_Dn
st_vec = fn_ww__setup__extract_st_vec_from_st_Dn( st_Dn );

% LSQEXP
st_v_shear_lsqexp_PF = struct;
st_v_shear_lsqexp_PF.v_U = st_res.fn__exp_fit( v_lsqexp_fit_param, st_Dn.v_zm );
st_v_shear_lsqexp_PF.v_dU = st_res.fn__exp_fit_d( v_lsqexp_fit_param, st_Dn.v_zm );
st_v_shear_lsqexp_PF.v_ddU = st_res.fn__exp_fit_dd( v_lsqexp_fit_param, st_Dn.v_zm );
[ st_r_shear_lsqexp_PF ] = fn_ww__setup__create_shear_r_st__vec( st_vec, st_v_shear_lsqexp_PF, st_p );

st_v_shear_lsqexp_APX = struct;
st_v_shear_lsqexp_APX.v_U = st_res.fn__exp_fit( v_lsqexp_fit_param, st_res.v_z_cc );
st_v_shear_lsqexp_APX.v_dU = st_res.fn__exp_fit_d( v_lsqexp_fit_param, st_res.v_z_cc );
st_v_shear_lsqexp_APX.v_ddU = st_res.fn__exp_fit_dd( v_lsqexp_fit_param, st_res.v_z_cc );
[ st_r_shear_lsqexp_APX ] = fn_ww__setup__create_shear_r_st__vec( st_vec_cc, st_v_shear_lsqexp_APX, st_p );

% OLS
st_v_shear_olscb_PF = struct;
st_v_shear_olscb_PF.v_U = st_fit_olscb_PF.v_zs_fit_U;
st_v_shear_olscb_PF.v_dU = st_fit_olscb_PF.v_zs_fit_dU;
st_v_shear_olscb_PF.v_ddU = st_fit_olscb_PF.v_zs_fit_ddU;
[ st_r_shear_olscb_PF ] = fn_ww__setup__create_shear_r_st__vec( st_vec, st_v_shear_olscb_PF, st_p );

st_v_shear_olscb_APX = struct;
st_v_shear_olscb_APX.v_U = st_fit_olscb_APX.v_zs_fit_U;
st_v_shear_olscb_APX.v_dU = st_fit_olscb_APX.v_zs_fit_dU;
st_v_shear_olscb_APX.v_ddU = st_fit_olscb_APX.v_zs_fit_ddU;
[ st_r_shear_olscb_APX ] = fn_ww__setup__create_shear_r_st__vec( st_vec_cc, st_v_shear_olscb_APX, st_p );

% TLS
st_v_shear_tlscb_PF = struct;
st_v_shear_tlscb_PF.v_U = st_fit_tlscb_PF.v_zs_fit_U;
st_v_shear_tlscb_PF.v_dU = st_fit_tlscb_PF.v_zs_fit_dU;
st_v_shear_tlscb_PF.v_ddU =  st_fit_tlscb_PF.v_zs_fit_ddU;
[ st_r_shear_tlscb_PF ] = fn_ww__setup__create_shear_r_st__vec( st_vec, st_v_shear_tlscb_PF, st_p );

st_v_shear_tlscb_APX = struct;
st_v_shear_tlscb_APX.v_U = st_fit_tlscb_APX.v_zs_fit_U;
st_v_shear_tlscb_APX.v_dU = st_fit_tlscb_APX.v_zs_fit_dU;
st_v_shear_tlscb_APX.v_ddU = st_fit_tlscb_APX.v_zs_fit_ddU;
[ st_r_shear_tlscb_APX ] = fn_ww__setup__create_shear_r_st__vec( st_vec_cc, st_v_shear_tlscb_APX, st_p );

% Bezier
st_v_shear_bezier_PF = struct;
st_v_shear_bezier_PF.v_U = st_fit_bezier_PF.v_y;
st_v_shear_bezier_PF.v_dU = st_fit_bezier_PF.v_dy;
st_v_shear_bezier_PF.v_ddU = st_fit_bezier_PF.v_ddy;
[ st_r_shear_bezier_PF ] = fn_ww__setup__create_shear_r_st__vec( st_vec, st_v_shear_bezier_PF, st_p );

st_v_shear_bezier_APX = struct;
st_v_shear_bezier_APX.v_U = st_fit_bezier_APX.v_y;
st_v_shear_bezier_APX.v_dU = st_fit_bezier_APX.v_dy;
st_v_shear_bezier_APX.v_ddU = st_fit_bezier_APX.v_ddy;
[ st_r_shear_bezier_APX ] = fn_ww__setup__create_shear_r_st__vec( st_vec_cc, st_v_shear_bezier_APX, st_p );


if ( b_do_extr )
    
    % OLS extr
    st_v_shear_olscb_PF_extr = struct;
    st_v_shear_olscb_PF_extr.v_U = st_fit_olscb_PF_extr.v_zs_fit_U;
    st_v_shear_olscb_PF_extr.v_dU = st_fit_olscb_PF_extr.v_zs_fit_dU;
    st_v_shear_olscb_PF_extr.v_ddU = st_fit_olscb_PF_extr.v_zs_fit_ddU;
    [ st_r_shear_olscb_PF_extr ] = fn_ww__setup__create_shear_r_st__vec( st_vec, st_v_shear_olscb_PF_extr, st_p );

    st_v_shear_olscb_APX_extr = struct;
    st_v_shear_olscb_APX_extr.v_U = st_fit_olscb_APX_extr.v_zs_fit_U;
    st_v_shear_olscb_APX_extr.v_dU = st_fit_olscb_APX_extr.v_zs_fit_dU;
    st_v_shear_olscb_APX_extr.v_ddU = st_fit_olscb_APX_extr.v_zs_fit_ddU;
    [ st_r_shear_olscb_APX_extr ] = fn_ww__setup__create_shear_r_st__vec( st_vec_cc, st_v_shear_olscb_APX_extr, st_p );
    
    % TLS extr
    st_v_shear_tlscb_PF_extr = struct;
    st_v_shear_tlscb_PF_extr.v_U = st_fit_tlscb_PF_extr.v_zs_fit_U;
    st_v_shear_tlscb_PF_extr.v_dU = st_fit_tlscb_PF_extr.v_zs_fit_dU;
    st_v_shear_tlscb_PF_extr.v_ddU = st_fit_tlscb_PF_extr.v_zs_fit_ddU;
    [ st_r_shear_tlscb_PF_extr ] = fn_ww__setup__create_shear_r_st__vec( st_vec, st_v_shear_tlscb_PF_extr, st_p );
    
    st_v_shear_tlscb_APX_extr = struct;
    st_v_shear_tlscb_APX_extr.v_U = st_fit_tlscb_APX_extr.v_zs_fit_U;
    st_v_shear_tlscb_APX_extr.v_dU = st_fit_tlscb_APX_extr.v_zs_fit_dU;
    st_v_shear_tlscb_APX_extr.v_ddU = st_fit_tlscb_APX_extr.v_zs_fit_ddU;
    [ st_r_shear_tlscb_APX_extr ] = fn_ww__setup__create_shear_r_st__vec( st_vec_cc, st_v_shear_tlscb_APX_extr, st_p );
    
    % OLS extr2
    st_v_shear_olscb_PF_extr2 = st_v_shear_olscb_PF_extr;
    st_v_shear_olscb_PF_extr2.v_U( v_hs_layer_idxs ) = v_lin_extr;
    st_v_shear_olscb_PF_extr2.v_dU( v_hs_layer_idxs ) = m_hs * ones( size( v_hs_layer_idxs ) );
    st_v_shear_olscb_PF_extr2.v_ddU( v_hs_layer_idxs ) = 0 * v_hs_layer_idxs;
    [ st_r_shear_olscb_PF_extr2 ] = fn_ww__setup__create_shear_r_st__vec( st_vec, st_v_shear_olscb_PF_extr2, st_p );
    
end

% Making sure we don't do anything stupid by later using the unscalled
% quantities... not that I'd do that.
clear st_fit_olscb_PF;
clear st_fit_tlscb_PF;
clear st_fit_bezier_PF;
clear st_fit_olscb_APX;
clear st_fit_tlscb_APX;
clear st_fit_bezier_APX;
if ( b_do_extr )
    clear st_fit_olscb_PF_extr;
    clear st_fit_tlscb_PF_extr;
    clear st_fit_olscb_APX_extr;
    clear st_fit_tlscb_APX_extr;
end


% Exp settings
%fn_local__do_simple_plot( 'figures_shear_prtb\fig__intro_profile_exp',  0, 1.0 );
% Powerlaw settings
%fn_local__do_simple_plot( 'figures_shear_prtb\fig__intro_profile_powerlaw', 0, 1.0 );
% Colombia settings
% fn_local__do_simple_plot( 'figures_shear_prtb\fig__intro_profile_columbia', -0.5, 1.0 );
% error()

% fn_local__do_fit_plots( st_v_shear_olscb_PF, 'figures_shear_prtb\fig__fit_exp_U__olscb_shear', 'figures_shear_prtb\fig__fit_exp_dU__olscb_shear', 'figures_shear_prtb\fig__fit_exp_ddU__olscb_shear', { 'a', 'b', 'c' } );
% fn_local__do_fit_plots( st_v_shear_tlscb_PF, 'figures_shear_prtb\fig__fit_exp_U__tlscb_shear', 'figures_shear_prtb\fig__fit_exp_dU__tlscb_shear', 'figures_shear_prtb\fig__fit_exp_ddU__tlscb_shear', { 'd', 'e', 'f' } );
% fn_local__do_fit_plots( st_v_shear_lsqexp_PF, 'figures_shear_prtb\fig__fit_exp_U__exp_shear', 'figures_shear_prtb\fig__fit_exp_dU__exp_shear', 'figures_shear_prtb\fig__fit_exp_ddU__exp_shear', { 'g', 'h', 'i' } );
% fn_local__do_fit_plots( st_v_shear_olscb_PF_extr, 'figures_shear_prtb\fig__fit_exp_U__olscb_shear_extr', 'figures_shear_prtb\fig__fit_exp_dU__olscb_shear_extr', 'figures_shear_prtb\fig__fit_exp_ddU__olscb_shear_extr', { 'j', 'k', 'l' } );
% %fn_local__do_fit_plots( st_v_shear_tlscb_PF_extr, 'figures_shear_prtb\fig__fit_exp_U__tlscb_shear_extr', 'figures_shear_prtb\fig__fit_exp_dU__tlscb_shear_extr', 'figures_shear_prtb\fig__fit_exp_ddU__tlscb_shear_extr', { 'm', 'n', 'o' } );
% error()
% 

% 
% figure(1);
% plot( st_Dn.v_zm, st_res.st_fn_shear_exact.fn_U( st_Dn.v_zm ), 'b' );
% hold on;
% plot( st_Dn.v_zm, st_v_shear_olscb_PF.v_U, 'r--' );
% plot( st_vec_cc.v_zm, st_v_shear_olscb_APX.v_U, 'm-.' );
% plot( v_zs, v_zs_candidate_U, 'r*' );
% hold off;
% 
% 
% figure(10);
% %plot( st_Dn.v_zm, st_res.st_fn_shear_exact.fn_U( st_Dn.v_zm ), 'b' );
% plot( st_Dn.v_zm, st_v_shear_olscb_PF.v_U - st_res.st_fn_shear_exact.fn_U( st_Dn.v_zm ), 'r--' );
% hold on;
% plot( st_vec_cc.v_zm, st_v_shear_olscb_APX.v_U - st_res.st_fn_shear_exact.fn_U( st_vec_cc.v_zm ), 'm-.' );
% plot( v_zs, v_zs_candidate_U - st_res.st_fn_shear_exact.fn_U( v_zs ), 'r*' );
% hold off;
% 
% 
% figure(2);
% plot( st_Dn.v_zm, st_res.st_fn_shear_exact.fn_dU( st_Dn.v_zm ), 'b' );
% hold on;
% plot( st_Dn.v_zm, st_v_shear_olscb_PF.v_dU, 'r--' );
% plot( st_vec_cc.v_zm, st_v_shear_olscb_APX.v_dU, 'm-.' );
% plot( st_Dn.v_zm, st_Dn.a_Dm * st_res.st_fn_shear_exact.fn_U( st_Dn.v_zm ), 'g-.', 'LineWidth', 1.5 );
% hold off;
% 
% 
% figure(3);
% plot( st_Dn.v_zm, st_res.st_fn_shear_exact.fn_ddU( st_Dn.v_zm ), 'b' );
% hold on;
% plot( st_Dn.v_zm, st_v_shear_olscb_PF.v_ddU, 'r--' );
% plot( st_vec_cc.v_zm, st_v_shear_olscb_APX.v_ddU, 'm-.' );
% plot( st_Dn.v_zm, st_Dn.a_D2m * st_res.st_fn_shear_exact.fn_U( st_Dn.v_zm ), 'g-.', 'LineWidth', 1.5 );
% hold off;
% 
% st_p
% 
% error();
% 





%--------------------------------------------------------------------------
% Solve dispersion relation using fitted shear profiles
%--------------------------------------------------------------------------

% We do the APX method first, as this allows us to test for critical layers
% and avoid a more expensive PF calculation when it won't work.

% Use APX method
[ v_c_lsqexp_APX ] = fn_ww__calc_re__apx_el_cc__red_c( v_k, st_res.v_z_cc, st_res.v_zw_cc, st_r_shear_lsqexp_APX, st_p );
[ v_c_olscb_APX ] = fn_ww__calc_re__apx_el_cc__red_c( v_k, st_res.v_z_cc, st_res.v_zw_cc, st_r_shear_olscb_APX, st_p );
[ v_c_tlscb_APX ] = fn_ww__calc_re__apx_el_cc__red_c( v_k, st_res.v_z_cc, st_res.v_zw_cc, st_r_shear_tlscb_APX, st_p );
[ v_c_bezier_APX ] = fn_ww__calc_re__apx_el_cc__red_c( v_k, st_res.v_z_cc, st_res.v_zw_cc, st_r_shear_bezier_APX, st_p );
if ( b_do_extr )
    [ v_c_olscb_APX_extr ] = fn_ww__calc_re__apx_el_cc__red_c( v_k, st_res.v_z_cc, st_res.v_zw_cc, st_r_shear_olscb_APX_extr, st_p );
    [ v_c_tlscb_APX_extr ] = fn_ww__calc_re__apx_el_cc__red_c( v_k, st_res.v_z_cc, st_res.v_zw_cc, st_r_shear_tlscb_APX_extr, st_p );
end


% Crit layer tests
b_crit_layer = false;
if ( min( v_c_lsqexp_APX ) < st_r_shear_lsqexp_APX.st_crit.U_max )
    b_crit_layer = true;
    if ( b_show_warnings ), warning( 'Crit layer in lsqexp APX' ); end
end
if ( min( v_c_olscb_APX ) < st_r_shear_olscb_APX.st_crit.U_max )
    b_crit_layer = true;
    if ( b_show_warnings ), warning( 'Crit layer in OLSCB APX' ); end
end
if ( min( v_c_tlscb_APX ) < st_r_shear_tlscb_APX.st_crit.U_max )
    b_crit_layer = true;
    if ( b_show_warnings ), warning( 'Crit layer in TLSCB APX' ); end
end
if ( min( v_c_bezier_APX ) < st_r_shear_bezier_APX.st_crit.U_max )
    b_crit_layer = true;
    if ( b_show_warnings ), warning( 'Crit layer in Bezier APX' ); end
end
if ( b_do_extr )
    if ( min( v_c_olscb_APX_extr ) < st_r_shear_olscb_APX_extr.st_crit.U_max )
        b_crit_layer = true;
        if ( b_show_warnings ), warning( 'Crit layer in OLSCB_extr APX' ); end
    end
    if ( min( v_c_tlscb_APX_extr ) < st_r_shear_tlscb_APX_extr.st_crit.U_max )
        b_crit_layer = true;
        if ( b_show_warnings ), warning( 'Crit layer in TLSCB_extr APX' ); end
    end        
end


% Set flag for later
if ( b_crit_layer )
    b_safe_result = false;
else
    b_safe_result = true;
end



if ( ~b_crit_layer ) 
    
    % Use ADPTV method
    b_pf_fail = false;
    st_proc = struct;
    st_proc.ip_method = 2;
    st_proc.pf_tol = st_res.pf_tol;
    try

        [ v_c_lsqexp_PF ] = fn_ww__calc_re__adptv__red_c( st_Dn, st_r_shear_lsqexp_PF, st_res.st_adptv, st_proc, st_p );
        [ v_c_olscb_PF ] = fn_ww__calc_re__adptv__red_c( st_Dn, st_r_shear_olscb_PF, st_res.st_adptv, st_proc, st_p );
        [ v_c_tlscb_PF ] = fn_ww__calc_re__adptv__red_c( st_Dn, st_r_shear_tlscb_PF, st_res.st_adptv, st_proc, st_p );
        [ v_c_bezier_PF ] = fn_ww__calc_re__adptv__red_c( st_Dn, st_r_shear_bezier_PF, st_res.st_adptv, st_proc, st_p );
        if ( b_do_extr )
            [ v_c_olscb_PF_extr ] = fn_ww__calc_re__adptv__red_c( st_Dn, st_r_shear_olscb_PF_extr, st_res.st_adptv, st_proc, st_p );
            [ v_c_tlscb_PF_extr ] = fn_ww__calc_re__adptv__red_c( st_Dn, st_r_shear_tlscb_PF_extr, st_res.st_adptv, st_proc, st_p );
            [ v_c_olscb_PF_extr2 ] = fn_ww__calc_re__adptv__red_c( st_Dn, st_r_shear_olscb_PF_extr2, st_res.st_adptv, st_proc, st_p );
        end

    catch obj_ex

        b_pf_fail = true;

    end

    %norm( v_c_bezier_PF - v_c_bezier_APX, inf )

    % plot( v_k, v_c_olscb_PF, 'b', v_k, v_c_olscb_APX, 'b:' );
    % hold on;
    % plot( v_k, v_c_bezier_PF, 'm', v_k, v_c_bezier_APX, 'm:' );
    % plot( v_k, v_c_lsqexp_PF, 'r', v_k, v_c_lsqexp_APX, 'r:' );
    % hold off;
    % error();



    % figure(1)
    % plot( v_k, v_c_lsqexp_PF, 'b', v_k, v_c_lsqexp_APX, 'b--' );
    % hold on;
    % plot( v_k, v_c_olscb_PF, 'm', v_k, v_c_olscb_APX, 'm--' );
    % plot( v_k, v_c_tlscb_PF, 'g', v_k, v_c_tlscb_APX, 'g--' );
    % hold off;
    % 
    % 
    % 
    % figure(10);
    % norm( v_c_lsqexp_PF - v_c_lsqexp_APX, inf )
    % norm( v_c_olscb_PF - v_c_olscb_APX, inf )
    % norm( v_c_tlscb_PF - v_c_tlscb_APX, inf )
    % norm( v_c_olscb_PF_extr - v_c_olscb_APX_extr, inf )
    % norm( v_c_tlscb_PF_extr - v_c_tlscb_APX_extr, inf )
    % plot( v_k, v_c_lsqexp_PF - v_c_lsqexp_APX, 'b' );
    % hold on;
    % plot( v_k, v_c_olscb_PF - v_c_olscb_APX, 'b' );
    % plot( v_k, v_c_tlscb_PF - v_c_tlscb_APX, 'b' );
    % plot( v_k, v_c_olscb_PF_extr - v_c_olscb_APX_extr, 'b' );
    % plot( v_k, v_c_tlscb_PF_extr - v_c_tlscb_APX_extr, 'b' );
    % error();



    % If the PF calculation went ok then run additional error checks
    if ( b_pf_fail )

        % We know it's bad
        b_safe_result = false;

    else

        % Set default
        b_safe_result = true;    

        % Check is monotonically decreasing
        if ( ~all( diff( v_c_lsqexp_PF ) < 0 ) )
            if ( b_show_warnings ), warning( 'LSQExp result not monotonically decreasing' ); end
            b_safe_result = false;
        end
        if ( ~all( diff( v_c_olscb_PF ) < 0 ) )
            if ( b_show_warnings ), warning( 'OLSCB result not monotonically decreasing' ); end
            b_safe_result = false;
        end
        if ( ~all( diff( v_c_tlscb_PF ) < 0 ) )
            if ( b_show_warnings ), warning( 'TLSCB result not monotonically decreasing' ); end
            b_safe_result = false;
        end
        % We don't use sanity tests on Bezier as it's too disruptive
%         if ( ~all( diff( v_c_bezier_PF ) < 0 ) )
%             if ( b_show_warnings ), warning( 'Bezier result not monotonically decreasing' ); end
%             b_safe_result = false;
%         end
        if ( b_do_extr )
            if ( ~all( diff( v_c_olscb_PF_extr ) < 0 ) )
                if ( b_show_warnings ), warning( 'OLSCB_extr result not monotonically decreasing' ); end
                b_safe_result = false;
            end
            if ( ~all( diff( v_c_tlscb_PF_extr ) < 0 ) )
                if ( b_show_warnings ), warning( 'TLSCB_extr result not monotonically decreasing' ); end
                b_safe_result = false;
            end
            if ( ~all( diff( v_c_olscb_PF_extr2 ) < 0 ) )
                if ( b_show_warnings ), warning( 'OLSCB_extr2 result not monotonically decreasing' ); end
                b_safe_result = false;
            end
        end

        % Check we're below a sane threshold
        if ( max( abs( ( v_c_lsqexp_PF - st_res.v_c_ref ) ./ st_res.v_c_ref ) ) > upperlim_relerr )
            if ( b_show_warnings ), warning( 'LSQExp result too large' ); end
            b_safe_result = false;
        end
        if ( max( abs( ( v_c_olscb_PF - st_res.v_c_ref ) ./ st_res.v_c_ref ) ) > upperlim_relerr )
            if ( b_show_warnings ), warning( 'OLSCB result too large' ); end
            b_safe_result = false;
        end
        if ( max( abs( ( v_c_tlscb_PF - st_res.v_c_ref ) ./ st_res.v_c_ref ) ) > upperlim_relerr )
            if ( b_show_warnings ), warning( 'TLSCB result too large' ); end
            b_safe_result = false;
        end
        % We don't use sanity tests on Bezier and it's too disruptive 
%         if ( max( abs( ( v_c_bezier_PF - st_res.v_c_ref ) ./ st_res.v_c_ref ) ) > upperlim_relerr )
%             if ( b_show_warnings ), warning( 'Bezier result too large' ); end
%             b_safe_result = false;
%         end
        if ( b_do_extr )
            if ( max( abs( ( v_c_olscb_PF_extr - st_res.v_c_ref ) ./ st_res.v_c_ref ) ) > upperlim_relerr )
                if ( b_show_warnings ), warning( 'OLSCB_extr result too large' ); end
                b_safe_result = false;
            end
            if ( max( abs( ( v_c_tlscb_PF_extr - st_res.v_c_ref ) ./ st_res.v_c_ref ) ) > upperlim_relerr )
                if ( b_show_warnings ), warning( 'TLSCB_extr result too large' ); end
                b_safe_result = false;
            end
            if ( max( abs( ( v_c_olscb_PF_extr2 - st_res.v_c_ref ) ./ st_res.v_c_ref ) ) > upperlim_relerr )
                if ( b_show_warnings ), warning( 'OLSCB_extr2 result too large' ); end
                b_safe_result = false;
            end
        end

        % Check PF isn't too different from APX
        if ( norm( v_c_lsqexp_PF - v_c_lsqexp_APX, inf ) > upperlim_clpf_apx_diff )
            if ( b_show_warnings ), warning( 'LSQExp PF and APX results too different' ); end
            b_safe_result = false;
        end
        if ( norm( v_c_olscb_PF - v_c_olscb_APX, inf ) > upperlim_clpf_apx_diff )
            if ( b_show_warnings ), warning( 'OLSCB PF and APX results too different' ); end
            b_safe_result = false;
        end
        if ( norm( v_c_tlscb_PF - v_c_tlscb_APX, inf ) > upperlim_clpf_apx_diff )
            if ( b_show_warnings ), warning( 'TLSCB PF and APX results too different' ); end
            b_safe_result = false;
        end
        % We don't use sanity tests on Bezier and it's too disruptive
%         if ( norm( v_c_bezier_PF - v_c_bezier_APX, inf ) > upperlim_clpf_apx_diff )
%             if ( b_show_warnings ), warning( 'Bezier PF and APX results too different' ); end
%             b_safe_result = false;
%         end
        if ( b_do_extr )
            if ( norm( v_c_olscb_PF_extr - v_c_olscb_APX_extr, inf ) > upperlim_clpf_apx_diff )
                if ( b_show_warnings ), warning( 'OLSCB_extr PF and APX results too different' ); end
                b_safe_result = false;
            end
            if ( norm( v_c_tlscb_PF_extr - v_c_tlscb_APX_extr, inf ) > upperlim_clpf_apx_diff )
                if ( b_show_warnings ), warning( 'TLSCB_extr PF and APX results too different' ); end
                b_safe_result = false;
            end
        end

    end    

end % End of PF initial, initial test of crit layer from APX





% Store results
st_prtb_single = struct;
if ( ~b_safe_result )
    
    st_err.b_ok = false;
%    warning( 'Hit unsafe result, only saving APX results shall zero PF' );
    
    st_prtb_single.v_c_lsqexp_PF = 0 * v_k;
    st_prtb_single.v_c_olscb_PF = 0 * v_k;
    st_prtb_single.v_c_tlscb_PF = 0 * v_k;
    st_prtb_single.v_c_bezier_PF = 0 * v_k;
    if ( b_do_extr )
        st_prtb_single.v_c_olscb_PF_extr = 0 * v_k;
        st_prtb_single.v_c_tlscb_PF_extr = 0 * v_k;
        st_prtb_single.v_c_olscb_PF_extr2 = 0 * v_k;
    end
    
else 

    st_prtb_single.v_c_lsqexp_PF = v_c_lsqexp_PF;
    st_prtb_single.v_c_olscb_PF = v_c_olscb_PF;
    st_prtb_single.v_c_tlscb_PF = v_c_tlscb_PF;
    st_prtb_single.v_c_bezier_PF = v_c_bezier_PF;
    if ( b_do_extr )
        st_prtb_single.v_c_olscb_PF_extr = v_c_olscb_PF_extr;
        st_prtb_single.v_c_tlscb_PF_extr = v_c_tlscb_PF_extr;
        st_prtb_single.v_c_olscb_PF_extr2 = v_c_olscb_PF_extr2;
    end

end


% Now save the APX results (we always save APX)
st_prtb_single.v_c_lsqexp_APX = v_c_lsqexp_APX;
st_prtb_single.v_c_olscb_APX = v_c_olscb_APX;
st_prtb_single.v_c_tlscb_APX = v_c_tlscb_APX;
st_prtb_single.v_c_bezier_APX = v_c_bezier_APX;
if ( b_do_extr )
    st_prtb_single.v_c_olscb_APX_extr = v_c_olscb_APX_extr;
    st_prtb_single.v_c_tlscb_APX_extr = v_c_tlscb_APX_extr;
end











   function fn_local__do_simple_plot( s_filename_simple, y_min, y_max )
        
%         y_min = -0.5;
%         y_max = 1.5;
        
        % Fontsizes
        fontsz_labels = 36;
        fontsz_axes = 32;

        % Linewidths
        line_width_std = 3.0;
        
        % Colours
        col_blue = [ 34 98 247 ] / 255; 
        col_blue_light = [ 160 184 239  ] / 255;

        col_green = [ 42 94 70 ] / 255;
        col_green_light = [ 184 224 206 ] / 255; 

        col_magenta = [ 173 12 214 ] / 255;
        col_magenta_light = [ 223 197 229 ] / 255;

        col_red = [ 239 52 52 ] / 255;
        col_red_light = [ 239 160 160 ] / 255;
        
        col_nearblack = [ 0.1 0.1 0.1 ];                
        
        
        % Simple U (no fitting or anything like that)
        h_fig_simple = figure(1);        
        plot( st_res.st_fn_shear_exact.fn_U( st_Dn.v_zm ), st_Dn.v_zm, 'Color', col_blue, 'LineWidth', line_width_std );
        ylim( [ -st_p.h 0 ] );
        xlim( [ y_min y_max ] );    
        %yticks( -st_p.h:2:0 );
        yticks( -st_p.h:0.25:0 );
        set(gca,'fontsize', fontsz_axes );
        set(gca,'TickLabelInterpreter', 'latex');        
        ylabel( '$\tilde{z}$', 'FontSize', fontsz_labels );
        xlabel( '$\tilde{U}(\tilde{z})$', 'FontSize', fontsz_labels );       
        hold off;
        
        set(gcf,'PaperOrientation','portrait');  %Changed
        set(gcf,'PaperSize', [ 20 22 ] );        
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPositionMode','auto');
        set(gcf,'PaperUnits','normalized');
        set(gcf,'PaperPosition', [0 0 1 1]);
        
        set( h_fig_simple, 'Color', 'none' );
        set( h_fig_simple, 'InvertHardCopy', 'Off' );        
        print( s_filename_simple, '-dpdf' );        

   end





    function fn_local__do_fit_plots( st_v_shear_plt, s_filename_U, s_filename_dU, s_filename_ddU, ca_lbls )
        
        y_min = -0.5;
        y_max = 1.5;
        
        % Fontsizes
        fontsz_labels = 38;
        fontsz_axes = 36;

        % Linewidths
        line_width_std = 3.0;
        
        % Colours
        col_blue = [ 34 98 247 ] / 255; 
        col_blue_light = [ 160 184 239  ] / 255;

        col_green = [ 42 94 70 ] / 255;
        col_green_light = [ 184 224 206 ] / 255; 

        col_magenta = [ 173 12 214 ] / 255;
        col_magenta_light = [ 223 197 229 ] / 255;

        col_red = [ 239 52 52 ] / 255;
        col_red_light = [ 239 160 160 ] / 255;
        
        col_nearblack = [ 0.1 0.1 0.1 ];
        
            
        % U
        h_fig1 = figure(1);        
        for lp_l=1:numel( v_zs )
            line( [ 1.25 1.75 ], [ v_zs(lp_l) v_zs(lp_l) ], 'LineStyle', ':', 'Color', 0.6 * ones( 1, 3 ), 'LineWidth', 1.5 );
            hold on;
        end
        line( [ 1.0 1.85 ], [ -hs -hs ], 'LineStyle', ':', 'Color', 0.3 * ones( 1, 3 ), 'LineWidth', 2.0 );
        line( [ y_min y_max ], [ 0 0 ], 'LineStyle', ':', 'Color', 0.6 * ones( 1, 3 ), 'LineWidth', 1.0 );
        plot( st_res.st_fn_shear_exact.fn_U( st_Dn.v_zm ), st_Dn.v_zm, 'Color', col_blue_light, 'LineWidth', line_width_std );        
        plot( st_v_shear_plt.v_U, st_Dn.v_zm, '--', 'Color', col_red, 'LineWidth', line_width_std );        
        plot( v_zs_candidate_U, v_zs, '*', 'Color', col_red, 'MarkerSize', 15, 'LineWidth', 2.0 );
        ylim( [ -st_p.h 0 ] );
        xlim( [ y_min y_max ] );    
        yticks( -st_p.h:2:0 );        
        set(gca,'fontsize', fontsz_axes );
        set(gca,'TickLabelInterpreter', 'latex');
        ylabel( '$\tilde{z}$', 'FontSize', fontsz_labels );
        xlabel( '$\tilde{U}(\tilde{z})$', 'FontSize', fontsz_labels );       
        hold off;       
        
        ob_title_a = title( sprintf( '(%s)', ca_lbls{1} ) );
        cur_title_pos = ob_title_a.Position
        cur_title_pos(1) = cur_title_pos(1) - 1.3;
        cur_title_pos(2) = cur_title_pos(2) - 0.5;
        set( ob_title_a, 'Position', cur_title_pos );

        
        set(gcf,'PaperOrientation','portrait');  %Changed
        set(gcf,'PaperSize', [ 20 22 ] );        
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPositionMode','auto');
        set(gcf,'PaperUnits','normalized');
        set(gcf,'PaperPosition', [0 0 1 1]);
        
        set( h_fig1, 'Color', 'none' );
        set( h_fig1, 'InvertHardCopy', 'Off' );        
        print( s_filename_U, '-dpdf' );
                

        % dU
        h_fig2 = figure(2);        
        for lp_l=1:numel( v_zs )
            line( [ 1.25 1.75 ], [ v_zs(lp_l) v_zs(lp_l) ], 'LineStyle', ':', 'Color', 0.6 * ones( 1, 3 ), 'LineWidth', 1.5 );
            hold on;
        end        
        line( [ 1.0 1.85 ], [ -hs -hs ], 'LineStyle', ':', 'Color', 0.3 * ones( 1, 3 ), 'LineWidth', 2.0 );
        line( [ y_min y_max ], [ 0 0 ], 'LineStyle', ':', 'Color', 0.6 * ones( 1, 3 ), 'LineWidth', 1.0 );
        plot( st_res.st_fn_shear_exact.fn_dU( st_Dn.v_zm ), st_Dn.v_zm, 'Color', col_blue_light, 'LineWidth', line_width_std );
        plot( st_v_shear_plt.v_dU, st_Dn.v_zm, '--', 'Color', col_red, 'LineWidth', line_width_std );
        ylim( [ -st_p.h 0 ] );
        xlim( [ y_min y_max ] );
        yticks( -st_p.h:2:0 );
        set(gca,'fontsize', fontsz_axes );
        set(gca,'TickLabelInterpreter', 'latex');
        ylabel( '$\tilde{z}$', 'FontSize', fontsz_labels );
        xlabel( '$\tilde{U}''(\tilde{z})$', 'FontSize', fontsz_labels );
        hold off;
        
        ob_title_b = title( sprintf( '(%s)', ca_lbls{2} ) );
        cur_title_pos = ob_title_b.Position
        cur_title_pos(1) = cur_title_pos(1) - 1.3;
        cur_title_pos(2) = cur_title_pos(2) - 0.5;
        set( ob_title_b, 'Position', cur_title_pos );
        
        set(gcf,'PaperOrientation','portrait');  %Changed
        set(gcf,'PaperSize', [ 20 22 ] );        
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPositionMode','auto');
        set(gcf,'PaperUnits','normalized');
        set(gcf,'PaperPosition', [0 0 1 1]);
        
        set( h_fig2, 'Color', 'none' );
        set( h_fig2, 'InvertHardCopy', 'Off' );
        print( s_filename_dU, '-dpdf' );

        
        % ddU
        h_fig3 = figure(3);
        for lp_l=1:numel( v_zs )
            line( [ 1.25 1.75 ], [ v_zs(lp_l) v_zs(lp_l) ], 'LineStyle', ':', 'Color', 0.6 * ones( 1, 3 ), 'LineWidth', 1.5 );
            hold on;
        end
        line( [ 1.0 1.85 ], [ -hs -hs ], 'LineStyle', ':', 'Color', 0.3 * ones( 1, 3 ), 'LineWidth', 2.0 );
        line( [ y_min y_max ], [ 0 0 ], 'LineStyle', ':', 'Color', 0.6 * ones( 1, 3 ), 'LineWidth', 1.0 );
        plot( st_res.st_fn_shear_exact.fn_ddU( st_Dn.v_zm ), st_Dn.v_zm, 'Color', col_blue_light, 'LineWidth', line_width_std );
        plot( st_v_shear_plt.v_ddU, st_Dn.v_zm, '--', 'Color', col_red, 'LineWidth', line_width_std );
        ylim( [ -st_p.h 0 ] );
        xlim( [ y_min y_max ] );
        yticks( -st_p.h:2:0 );
        set(gca,'fontsize', fontsz_axes );
        set(gca,'TickLabelInterpreter', 'latex');
        ylabel( '$\tilde{z}$', 'FontSize', fontsz_labels );
        xlabel( '$\tilde{U}''''(\tilde{z})$', 'FontSize', fontsz_labels );
        hold off;
        
        ob_title_c = title( sprintf( '(%s)', ca_lbls{3} ) );
        cur_title_pos = ob_title_c.Position
        cur_title_pos(1) = cur_title_pos(1) - 1.3;
        cur_title_pos(2) = cur_title_pos(2) - 0.5;
        set( ob_title_c, 'Position', cur_title_pos );
        
        set(gcf,'PaperOrientation','portrait');  %Changed
        set(gcf,'PaperSize', [ 20 22 ] );
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPositionMode','auto');
        set(gcf,'PaperUnits','normalized');
        set(gcf,'PaperPosition', [0 0 1 1]);
        
        set(h_fig3, 'Color', 'none' );
        set(h_fig3, 'InvertHardCopy', 'Off' );        
        print( s_filename_ddU, '-dpdf' );

        
        pause(2);
        close( h_fig1 );
        close( h_fig2 );
        close( h_fig3 );

    end





                








% 
%     function fn_local__do_fit_plots( st_v_shear_plt, s_filename_U, s_filename_dU, s_filename_ddU )
%         
%         y_min = -0.5;
%         y_max = 2.0;
%         
%         % Fontsizes
%         fontsz_labels = 26;
%         fontsz_axes = 22;
% 
%         % Linewidths
%         line_width_std = 3.0;
%         
%         % Colours
%         col_blue = [ 34 98 247 ] / 255; 
%         col_blue_light = [ 160 184 239  ] / 255;
% 
%         col_green = [ 42 94 70 ] / 255;
%         col_green_light = [ 184 224 206 ] / 255; 
% 
%         col_magenta = [ 173 12 214 ] / 255;
%         col_magenta_light = [ 223 197 229 ] / 255;
% 
%         col_red = [ 239 52 52 ] / 255;
%         col_red_light = [ 239 160 160 ] / 255;
%         
%         col_nearblack = [ 0.1 0.1 0.1 ];
%         
%         % U
%         h_fig1 = figure(1);        
%         for lp_l=1:numel( v_zs )
%             line( [ v_zs(lp_l) v_zs(lp_l) ],  [ 1.25 1.75 ], 'LineStyle', ':', 'Color', 0.6 * ones( 1, 3 ), 'LineWidth', 2.0 );
%             hold on;
%         end
%         line( [ -hs -hs ], [ 1.15 1.85 ], 'LineStyle', ':', 'Color', 0.3 * ones( 1, 3 ), 'LineWidth', 2.5 );
%         plot( st_Dn.v_zm, st_res.st_fn_shear_exact.fn_U( st_Dn.v_zm ), 'Color', col_blue_light, 'LineWidth', line_width_std );        
%         plot( st_Dn.v_zm, st_v_shear_plt.v_U , '--', 'Color', col_red, 'LineWidth', line_width_std );
%         %plot( v_zs, -0.5 + 0*v_zs, 'mx', 'MarkerSize', 10 );
%         plot( v_zs, v_zs_candidate_U, '*', 'Color', col_red, 'MarkerSize', 15, 'LineWidth', 2.0 );
%         xlabel( '$\tilde{z}$', 'FontSize', fontsz_labels );
%         ylabel( '$\tilde{U}(\tilde{z})$', 'FontSize', fontsz_labels );
%         xlim( [ -st_p.h 0 ] );
%         ylim( [ y_min y_max ] );    
%         xticks( -st_p.h:1:0 );
%         set(gca,'fontsize', fontsz_axes );
%         set(gca,'TickLabelInterpreter', 'latex');
%         hold off;
%         set(gcf,'PaperOrientation','landscape');
%         set(gcf,'PaperPositionMode','auto');
%         set(gcf,'PaperUnits','normalized');
%         set(gcf,'PaperPosition', [0 0 1 1]);
%         print( s_filename_U, '-dpdf' );
%             
%         
%         % dU
%         h_fig2 = figure(2);        
%         for lp_l=1:numel( v_zs )
%             line( [ v_zs(lp_l) v_zs(lp_l) ],  [ 1.25 1.75 ], 'LineStyle', ':', 'Color', 0.6 * ones( 1, 3 ), 'LineWidth', 2.0 );
%             hold on;
%         end        
%         line( [ -hs -hs ], [ 1.15 1.85 ], 'LineStyle', ':', 'Color', 0.3 * ones( 1, 3 ), 'LineWidth', 2.5 );
%         plot( st_Dn.v_zm, st_res.st_fn_shear_exact.fn_dU( st_Dn.v_zm ), 'Color', col_blue_light, 'LineWidth', line_width_std );
%         plot( st_Dn.v_zm, st_v_shear_plt.v_dU , '--', 'Color', col_red, 'LineWidth', line_width_std );
%         xlabel( '$\tilde{z}$', 'FontSize', fontsz_labels );
%         ylabel( '$\tilde{U}''(\tilde{z})$', 'FontSize', fontsz_labels );
%         xlim( [ -st_p.h 0 ] );
%         ylim( [ y_min y_max ] );
%         xticks( -st_p.h:1:0 );
%         set(gca,'fontsize', fontsz_axes );
%         set(gca,'TickLabelInterpreter', 'latex');
%         hold off;
%         set(gcf,'PaperOrientation','landscape');
%         set(gcf,'PaperPositionMode','auto');
%         set(gcf,'PaperUnits','normalized');
%         set(gcf,'PaperPosition', [0 0 1 1]);
%         print( s_filename_dU, '-dpdf' );
% 
%         
%         % ddU
%         h_fig3 = figure(3);
%         for lp_l=1:numel( v_zs )
%             line( [ v_zs(lp_l) v_zs(lp_l) ],  [ 1.25 1.75 ], 'LineStyle', ':', 'Color', 0.6 * ones( 1, 3 ), 'LineWidth', 2.0 );
%             hold on;
%         end
%         line( [ -hs -hs ], [ 1.15 1.85 ], 'LineStyle', ':', 'Color', 0.3 * ones( 1, 3 ), 'LineWidth', 2.5 );
%         plot( st_Dn.v_zm, st_res.st_fn_shear_exact.fn_ddU( st_Dn.v_zm ), 'Color', col_blue_light, 'LineWidth', line_width_std );
%         plot( st_Dn.v_zm, st_v_shear_plt.v_ddU, '--', 'Color', col_red, 'LineWidth', line_width_std );
%         xlabel( '$\tilde{z}$', 'FontSize', fontsz_labels );
%         ylabel( '$\tilde{U}''''(\tilde{z})$', 'FontSize', fontsz_labels );
%         xlim( [ -st_p.h 0 ] );
%         ylim( [ y_min y_max ] );
%         xticks( -st_p.h:1:0 );
%         set(gca,'fontsize', fontsz_axes );
%         set(gca,'TickLabelInterpreter', 'latex');
%         hold off;
%         set(gcf,'PaperOrientation','landscape');
%         set(gcf,'PaperPositionMode','auto');
%         set(gcf,'PaperUnits','normalized');
%         set(gcf,'PaperPosition', [0 0 1 1]);
%         print( s_filename_ddU, '-dpdf' );
% 
%         
%         pause(2);
%         close( h_fig1 );
%         close( h_fig2 );
%         close( h_fig3 );
% 
%     end
% 



end
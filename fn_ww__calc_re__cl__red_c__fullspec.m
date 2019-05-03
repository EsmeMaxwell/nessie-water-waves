function [ a_c_CL_pn, a_c_CL_safe ] = fn_ww__calc_re__cl__red_c__fullspec( st_Dn, v_k, st_r_shear, st_p )
%fn_ww__calc_re__cl__red_c__fullspec: INCOMPLETE  Calc dispersion relation for c + full spectrum
%

warning( 'This is a bit keech' );

% Basic sanity checks
assert( isstruct( st_Dn ) );
assert( isnumeric( v_k ) && 1 == size( v_k, 1 ) );
assert( isstruct( st_r_shear ) );
assert( isstruct( st_p ) );

% Set required mp digits (remember, it doesn't work well in parfor loops)
fn_ww__setup__prepare_mp( st_Dn, st_p );



N = numel( st_Dn.v_z0 );
Nk = numel( v_k );


% Get shear data
[ st_v_shear, U_nd_min, U_nd_max ] = fn_ww__util__get_r_shear_data( st_Dn, st_r_shear, st_p );

a_U = diag( st_v_shear.v_U + st_p.fp_current_shift );
a_dU = diag( st_v_shear.v_dU );
a_ddU = diag( st_v_shear.v_ddU );

U_min = U_nd_min + st_p.fp_current_shift;
U_max = U_nd_max + st_p.fp_current_shift;

if ( st_p.bp_mp )
    a_I = eye( size( a_U ), 'mp' );
else
    a_I = eye( size( a_U ) );
end



% Precalc
[ a_A0_precomp, a_A2, a_A0_FS, a_A1_FS, a_A2_FS ] = fn_ww__calc_re__cl__mtrxs__precalc_c( st_Dn, a_U, a_dU, a_ddU, a_I, st_p );

% Allocate memory
if ( st_p.bp_mp )
    a_c_CL_pn = zeros( 2, Nk, 'mp' );
    a_c_CL_safe = zeros( 2*N, Nk, 'mp' );
    
else
    a_c_CL_pn = zeros( 2, Nk );
    a_c_CL_safe = zeros( 2*N, Nk );
end



% Main calculation loop over v_k
for lp_k=1:Nk           
    
    if ( st_p.bp_disp_update ), fprintf( '%d ', lp_k ); end
    if ( st_p.bp_disp_update && 0 == mod( lp_k, 20 ) ), fprintf( '\n' ); end
    
    k = v_k(lp_k);
    
	% Core equations        
    [ a_A0, a_A1, a_A2 ] = fn_ww__calc_re__cl__mtrxs__core_c( st_Dn, a_U, a_A0_precomp, a_A2, a_I, a_A0_FS, a_A1_FS, a_A2_FS, k, st_p );
    
    % Possibly do FLD rescaling
    if ( st_p.bp_fld )
        [ a_A2, a_A1, a_A0, alpha ] = fn_ww__util__fld_scaling( a_A2, a_A1, a_A0 );
    end
    
    % Do the eigenvalue calculation
    [ a_w, v_eig ] = polyeig( a_A0, a_A1, a_A2 );
       
    % May have to put the eigenvalues back into the correct scale
    if ( st_p.bp_fld )
        v_eig = v_eig * alpha;
    end

    % Filter nonsense eigenvalue out
    [ v_eig, a_w, v_valid_idxs ] = fn_ww__calc_re__filter_valid_eigenvalues( v_eig, a_w, st_p );
    
    % Sort
    [ v_eig, v_sort_idx ] = sort( v_eig );
    a_w = a_w( :, v_sort_idx );
    
    a_c_CL_safe( 1:numel(v_eig), lp_k ) = v_eig;    
    
    % Find min and max
    [ c_max, ~ ] = fn_ww__calc_re__find_valid_eigvecs( st_Dn.v_zm, v_eig, a_w, U_min, U_max, 0 );
    [ c_min, ~ ] = fn_ww__calc_re__find_valid_eigvecs( st_Dn.v_zm, v_eig, a_w, U_min, U_max, 1 );
    a_c_CL_pn( 1, lp_k ) = c_max;
    a_c_CL_pn( 2, lp_k ) = c_min;
    

end
if ( st_p.bp_disp_update ), fprintf( '\n' ); end


% 
%     
%     % NEW code to find our desired eigenvectors directly... oooh, fancy.
%     
%     if ( 1 == lp_k )        
%         % Initialise    
%         v_eig_filter_inf_to_nan = v_eig;
%         v_eig_filter_inf_to_nan( ~isfinite( v_eig_filter_inf_to_nan ) ) = NaN;
%         v_eig_filter_inf_to_nan( abs( v_eig_filter_inf_to_nan ) > 1e3 ) = NaN;
%         
%         [ ~, max_idx ] = max( real( v_eig_filter_inf_to_nan ), [], 'omitnan' );
%         [ ~, min_idx ] = min( real( v_eig_filter_inf_to_nan ), [], 'omitnan' );
%         v_w_last_p = a_w( :, max_idx );
%         v_w_last_n = a_w( :, min_idx );
%        
%         a_c_CL_pn( 1, lp_k ) = v_eig( max_idx );
%         a_c_CL_pn( 2, lp_k ) = v_eig( min_idx );        
%         
%     else
%         fin_max_idx = -1;
%         fin_min_idx = -1;
%         max_opt_norm = 1e50;
%         min_opt_norm = 1e50;        
%                
%         for lp_ev=1:size( a_w, 2 )
%             
%             test_norm_max = norm( a_w( :, lp_ev ) - v_w_last_p, 1 );
%             if ( test_norm_max < max_opt_norm && isfinite( v_eig( lp_ev ) ) && abs( v_eig( lp_ev ) ) < 1e3 && abs( v_eig( lp_ev ) - a_c_CL_pn( 1, lp_k-1 ) ) < 0.15 )
%                 fin_max_idx = lp_ev;
%                 max_opt_norm = test_norm_max;
%             end
%             
%             test_norm_min = norm( a_w( :, lp_ev ) - v_w_last_n, 1 );
%             if ( test_norm_min < min_opt_norm && isfinite( v_eig( lp_ev ) ) && abs( v_eig( lp_ev ) ) < 1e3 && abs( v_eig( lp_ev ) - a_c_CL_pn( 2, lp_k-1 ) ) < 0.15 )
%                 fin_min_idx = lp_ev;
%                 min_opt_norm = test_norm_min;
%             end
%                         
%         end
%                        
%         v_w_last_p = a_w( :, fin_max_idx );
%         v_w_last_n = a_w( :, fin_min_idx );
%         
%         a_c_CL_pn( 1, lp_k ) = v_eig( fin_max_idx );
%         a_c_CL_pn( 2, lp_k ) = v_eig( fin_min_idx );
%         
%     end
%         
%     
%     % Now remove eigenvectors corresponding to dunted eigenvalues
%     v_nonfinite_idxs = find( ~isfinite( v_eig ) ).';
%     v_unsafe_idxs = find( abs( v_eig ) > st_p.fp_safety_thrsh ).';
%     v_imag_idxs = find( abs( imag( v_eig ) ) > st_p.fp_imag_thrsh ).';
%     
%     v_bad_e_idxs = union( v_nonfinite_idxs, v_unsafe_idxs );
%     v_bad_e_idxs = union( v_bad_e_idxs, v_imag_idxs );
%     v_valid_e_idxs = 1:numel(v_eig);
%     v_valid_e_idxs = setdiff( v_valid_e_idxs, v_bad_e_idxs );
%       
%     
%     % Possibly do extra filtering using eigenvectors
%     if ( false )        
%         [ v_valid_ev_idxs, ~ ] = fn_ww__calc_re__find_valid_eigvecs( st_Dn.v_zm, a_w, v_eig );
%         
%         v_bad_ev_idxs = 1:numel(v_eig);
%         v_bad_ev_idxs = setdiff( v_bad_ev_idxs, v_valid_ev_idxs );     
%         v_bad_ev_idxs
%     else
%         v_bad_ev_idxs = [];
%     end
%     
%     v_eig_safe = v_eig( v_valid_e_idxs );
%     v_eig_unsafe_e = v_eig( v_bad_e_idxs );
%     v_eig_unsafe_ev = v_eig( v_bad_ev_idxs );
%     
%     a_c_CL_safe( 1:numel(v_eig_safe), lp_k ) = v_eig_safe;
%     a_c_CL_unsafe_e( 1:numel(v_eig_unsafe_e), lp_k ) = v_eig_unsafe_e;
%     a_c_CL_unsafe_ev( 1:numel(v_eig_unsafe_ev), lp_k ) = v_eig_unsafe_ev;
%         
%     a_c_CL_safe( :, lp_k ) = sort( a_c_CL_safe( :, lp_k ), 'descend', 'ComparisonMethod', 'real' );
%     a_c_CL_unsafe_e( :, lp_k ) = sort( a_c_CL_unsafe_e( :, lp_k ), 'descend', 'ComparisonMethod', 'real' );
%     a_c_CL_unsafe_ev( :, lp_k ) = sort( a_c_CL_unsafe_ev( :, lp_k ), 'descend', 'ComparisonMethod', 'real' );
%     
%     
    
%     % Possibly delete zero eigenvalue (and vectors)
%     if ( 0 ~= st_p.bp_delete_zero_eig )
%         v_zero_idxs = find( abs( v_eig ) < st_p.fp_zero_tol );
%         v_eig( v_zero_idxs ) = [];
%         a_w( :, v_zero_idxs ) = [];
%     end    
    
%     % Error check
%     assert( numel( v_eig ) > 0, 'No valid eigenvalues found after filtering.' );
    

end
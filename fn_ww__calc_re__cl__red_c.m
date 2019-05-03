function [ v_c_p, a_c_w_p, v_residual, v_berr, v_cond ] = fn_ww__calc_re__cl__red_c( st_Dn, v_k, st_r_shear, st_p )
%fn_ww__calc_re__cl__red_c: Calc CL-c disperion relation for reduced profile (solve Rayleigh eqn with free-surface)
% 
%   [ v_c_p, a_c_w_p, v_residual, v_berr, v_cond ] = fn_ww__calc_re__cl__red_c( st_Dn, v_k, st_v_shear, st_p )
% 
% Calculates the dispersion relation and eigenvectors. Uses k as parameter
% and c as eigenvalue.
% 
% INPUT
%   st_Dn : Differentiation matrices
%   v_k : Vector of wave numbers
%   st_r_shear : Struct of reduced shear info
%   st_p : Parameters
% 
% OUTPUT
%   v_c_p : Vector of c results,  positive branch
%   a_c_w_p : Corresponding eigenvectors, arranged columnwise in k
%   (optional)
%   v_residual : Vector of residuals (optional)
%   v_berr : Vector of backwards errors (optional) 
%   v_cond : Vector of condition numbers (optional)
% 
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__diffmtrx__WR_poldif(),
%   fn_ww__setup__lin_map_Dn_to_mapped(),
%   fn_ww__setup__diffmtrx_mp__WR_chebdif(),
%   fn_ww__setup__shear_fn_to_vec(),
%   fn_ww__util__create_k_vec(),
%   fn_ww__calc_re__filter_valid_eigenvalues(),
%   fn_ww__calc_re__find_valid_eigvecs()


% TODO must convert all other CL routines to exclude bottom boundary condition

% TODO check Fr2 nonzero somewhere
% TODO check in calcs that st_Dn is matched to st_p and shear profile st_v_shear
% TODO in manuscript put emphasis re collocation on conditioning and finding correct eigenvectors, post-processing, etc

% TODO find U_min, U_max should be done analytically as vector may not be good enough

% TODO backport multicore stuff to here fn_ww__dev_re__cl_red_c__parfor()


% Basic sanity checks
assert( isstruct( st_Dn ) );
assert( isnumeric( v_k ) && 1 == size( v_k, 1 ) );
assert( isstruct( st_r_shear ) );
assert( isstruct( st_p ) );

assert( st_p.fp_current_shift == 0, 'Need to fix wrt crit layer' );




% Set required mp digits (remember, it doesn't work well in parfor loops)
fn_ww__setup__prepare_mp( st_Dn, st_p );

N = numel( st_Dn.v_z0 );
Nk = numel( v_k );



% U matrices and check crit layer position
if ( st_p.bp_phy_calc )        
    
    % Get shear data
    [ st_v_shear, U_nd_min, U_nd_max, U_phy_min, U_phy_max ] = fn_ww__util__get_r_shear_data( st_Dn, st_r_shear, st_p );    
    
    a_U = diag( st_v_shear.v_phy_U + st_p.fp_current_shift );
    a_dU = diag( st_v_shear.v_phy_dU );
    a_ddU = diag( st_v_shear.v_phy_ddU );
    
    U_min = U_phy_min + st_p.fp_current_shift;
    U_max = U_phy_max + st_p.fp_current_shift;
    
else
    
    % Get shear data
    [ st_v_shear, U_nd_min, U_nd_max ] = fn_ww__util__get_r_shear_data( st_Dn, st_r_shear, st_p );
    
    a_U = diag( st_v_shear.v_U + st_p.fp_current_shift );
    a_dU = diag( st_v_shear.v_dU );
    a_ddU = diag( st_v_shear.v_ddU );
    
    U_min = U_nd_min + st_p.fp_current_shift;
    U_max = U_nd_max + st_p.fp_current_shift;
    
end
if ( st_p.bp_mp )
    a_I = eye( size( a_U ), 'mp' );
else
    a_I = eye( size( a_U ) );
end




% Precalc
[ a_A0_precomp, a_A2, a_A0_FS, a_A1_FS, a_A2_FS ] = fn_ww__calc_re__cl__mtrxs__precalc_c( st_Dn, a_U, a_dU, a_ddU, a_I, st_p );


% Allocate memory
if ( st_p.bp_mp )
    v_c_p = zeros( 1, Nk, 'mp' );
else
    v_c_p = zeros( 1, Nk );
end

if ( nargout > 1 )
    if ( st_p.bp_mp )
        a_c_w_p = zeros( N, Nk, 'mp' );        
    else
        a_c_w_p = zeros( N, Nk );
    end
end

if ( nargout > 2 )
    if ( st_p.bp_mp )
        v_residual = zeros( 1, Nk, 'mp' );
        v_berr = zeros( 1, Nk, 'mp' );
    else
        v_residual = zeros( 1, Nk );
        v_berr = zeros( 1, Nk );
    end
    
    % Darn we're calculating condition number
    if ( nargout > 4 )
        if ( st_p.bp_mp )
            v_cond = zeros( 1, Nk, 'mp' );
        else
            v_cond = zeros( 1, Nk );
        end
    end    
    
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
    if ( st_p.bp_ex_btm_bc )
        
        % Calculate on rows without bottom BC
        if ( nargout > 4 )                        
            [ a_w, a_wl, v_eig ] = fn_ww__dev__polyeig( a_A0(1:end-1,1:end-1), a_A1(1:end-1,1:end-1), a_A2(1:end-1,1:end-1) );
             a_wl = [ a_wl; zeros( 1, size( a_w, 2 ) ) ];
             a_w = [ a_w; zeros( 1, size( a_w, 2 ) ) ];
        else
        	[ a_w, v_eig ] = polyeig( a_A0(1:end-1,1:end-1), a_A1(1:end-1,1:end-1), a_A2(1:end-1,1:end-1) );
            a_w = [ a_w; zeros( 1, size( a_w, 2 ) ) ];
        end
        
    else
        
        if ( nargout > 4 )                        
            [ a_w, a_wl, v_eig ] = fn_ww__dev__polyeig( a_A0, a_A1, a_A2 );
        else
            [ a_w, v_eig ] = polyeig( a_A0, a_A1, a_A2 );
        end
        
    end
       
    % May have to put the eigenvalues back into the correct scale
    if ( st_p.bp_fld )
        v_eig = v_eig * alpha;
    end
    
    % Now remove eigenvectors corresponding to dunted eigenvalues
    [ v_eig, a_w, v_valid_idxs ] = fn_ww__calc_re__filter_valid_eigenvalues( v_eig, a_w, st_p );           
    
    % Error check
    assert( numel( v_eig ) > 0, 'No valid eigenvalues found after eigenvalue filtering.' );    

    % Check for critical layer or zero eigenvalue, and process accordingly
    [ max_eig, max_idx ] = max( v_eig );
    if ( max_eig <= U_max || abs( max_eig ) < st_p.fp_zero_tol )
    
        warning( 'Hit critical layer' );
        
        % We're in critical layer terriroty, must be smart
        [ v_c_p(lp_k), a_c_w_p(:,lp_k) ] = fn_ww__calc_re__find_valid_eigvecs( st_Dn.v_zm, v_eig, a_w, U_min, U_max );
        
    else 
   
        % We're dominant, so can just use the maximum
        v_c_p(lp_k) = max_eig;
        
        % Pull eigenvector and normalise (and make sure sign consistent)
        if( nargout > 1 )            
            a_c_w_p( :, lp_k ) = fn_ww__util__normalise_eigenvectors( a_w( :, max_idx ), st_p.ip_std_evec_norm );

            % If we need to calculate residual
            if ( nargout > 2 )
                
                % We need the 2-norm for this...
                v_e = fn_ww__util__normalise_eigenvectors( a_w( :, max_idx ), 2 );
                
                % Calc residual
                v_residual( lp_k ) = norm( ( a_A2 * v_c_p(lp_k)^2 + a_A1 * v_c_p(lp_k) + a_A0 ) * v_e, 2 );
                
                % Calc sum of product of matrix eig norms
                sum_eig_mtrx_norm = norm( a_A2, 2 ) * abs( v_c_p(lp_k) )^2 + norm( a_A1, 2 ) * abs( v_c_p(lp_k) ) + norm( a_A0, 2 );
                
                % Calc backwards error (Tisseur thm 1 or Higham, Tissuer,
                % 2007)
                v_berr( lp_k ) = v_residual( lp_k ) / ( sum_eig_mtrx_norm * norm( v_e, 2 ) );
                
                % If we also need to calculate cond num
                if ( nargout > 4 )
                    
                    % Left eigenvector
                    v_el = fn_ww__util__normalise_eigenvectors( a_wl( :, max_idx ), 2 );
                    
                    % P(c) derivative
                    a_Pd = 2 *  v_c_p(lp_k) * a_A2 + a_A1;
                    
                    % Calc condition number
                    v_cond( lp_k ) = sum_eig_mtrx_norm * norm( v_el, 2 ) * norm( v_e, 2 ) / ( abs( v_c_p(lp_k) ) * abs( v_el' * a_Pd * v_e ) );
                    
                end                
                
            end

        end        
        
    end
        
    % Undo current shift, if required    
    v_c_p(lp_k) = v_c_p(lp_k) - st_p.fp_current_shift;    
    
end
if ( st_p.bp_disp_update ), fprintf( '\n' ); end




end
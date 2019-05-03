function [ v_c_p, a_c_w_p ] = fn_ww__calc_re__cl__gen_sl_pol_ang_c( st_Dn, k, v_theta, st_v_shear_x, st_v_shear_y, st_p )
%fn_ww__calc_re__cl__gen_sl_pol_ang_c: Calc CL-c disperion relation for angular general problem
% 
%   [ v_c_p, a_c_w_p ] = fn_ww__calc_re__cl__gen_sl_pol_ang_c( st_Dn, k, v_theta, st_v_shear_x, st_v_shear_y, st_p )
% 
% Calculates the dispersion relation and eigenvectors. Uses theta as
% parameter, fixed k, and c as eigenvalue.
% 
% INPUT
%   st_Dn : Differentiation matrices
%   k : Fixed wave number (radius)
%   v_theta : Vector of wave numbers
%   st_v_shear_x, st_v_shear_y : Struct of vectorised shear profile
%   st_p : Parameters
% 
% OUTPUT
%   v_c_p : Vector of c results,  positive branch
%   a_c_w_p : Corresponding eigenvectors, arranged columnwise in theta
%   (optional)
% 
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__diffmtrx__WR_poldif(),
%   fn_ww__setup__lin_map_Dn_to_mapped(),
%   fn_ww__setup__diffmtrx_mp__WR_chebdif(),
%   fn_ww__setup__shear_fn_to_vec()




% Basic sanity checks
assert( isstruct( st_Dn ) );
assert( isnumeric( k ) && 1 == numel( k ) );
assert( isnumeric( v_theta ) && 1 == size( v_theta, 1 ) );
assert( isstruct( st_v_shear_x ) );
assert( isstruct( st_v_shear_y ) );
assert( isstruct( st_p ) );

% Set required mp digits (remember, it doesn't work well in parfor loops)
fn_ww__setup__prepare_mp( st_Dn, st_p );




N = numel( st_Dn.v_z0 );
Ntheta = numel( v_theta );


if ( st_p.bp_mp )
    a_I = eye( size( st_Dn.a_Dm ), 'mp' );
else
    a_I = eye( size( st_Dn.a_Dm ) );
end


% Allocate memory
if ( st_p.bp_mp )
    v_c_p = zeros( 1, Ntheta, 'mp' );
else
    v_c_p = zeros( 1, Ntheta );
end

if ( nargout > 1 )
    if ( st_p.bp_mp )
        a_c_w_p = zeros( N, Ntheta, 'mp' );        
    else
        a_c_w_p = zeros( N, Ntheta );
    end
end




% Main calculation loop over v_theta
for lp_theta=1:Ntheta
        
    if ( st_p.bp_disp_update ), fprintf( '%d ', lp_theta ); end
    if ( st_p.bp_disp_update && 0 == mod( lp_theta, 20 ) ), fprintf( '\n' ); end    
    
    theta = v_theta(lp_theta);

    % U matrices
    if ( st_p.bp_phy_calc )
        v_phy_U = cos(theta) * st_v_shear_x.v_phy_U + sin(theta) * st_v_shear_y.v_phy_U + st_p.fp_current_shift;
        v_phy_dU = cos(theta) * st_v_shear_x.v_phy_dU + sin(theta) * st_v_shear_y.v_phy_dU;
        v_phy_ddU = cos(theta) * st_v_shear_x.v_phy_ddU + sin(theta) * st_v_shear_y.v_phy_ddU;                
        
        a_U = diag( v_phy_U + st_p.fp_current_shift );
        a_dU = diag( v_phy_dU );
        a_ddU = diag( v_phy_ddU );
        
        U_min = st_v_shear_x.U_phy_min + st_p.fp_current_shift;
        U_max = st_v_shear_x.U_phy_max + st_p.fp_current_shift;
    else
        v_U = cos(theta) * st_v_shear_x.v_U + sin(theta) * st_v_shear_y.v_U + st_p.fp_current_shift;
        v_dU = cos(theta) * st_v_shear_x.v_dU + sin(theta) * st_v_shear_y.v_dU;
        v_ddU = cos(theta) * st_v_shear_x.v_ddU + sin(theta) * st_v_shear_y.v_ddU;
        
        a_U = diag( v_U + st_p.fp_current_shift );
        a_dU = diag( v_dU );
        a_ddU = diag( v_ddU );
        
        U_min = st_v_shear.U_min + st_p.fp_current_shift;
        U_max = st_v_shear.U_max + st_p.fp_current_shift;
    end
    
    % Precalc
    [ a_A0_precomp, a_A2, a_A0_FS, a_A1_FS, a_A2_FS ] = fn_ww__calc_re__cl__mtrxs__precalc_c( st_Dn, a_U, a_dU, a_ddU, a_I, st_p );
        
    % Core equations
    [ a_A0, a_A1, a_A2 ] = fn_ww__calc_re__cl__mtrxs__core_c( st_Dn, a_U, a_A0_precomp, a_A2, a_I, a_A0_FS, a_A1_FS, a_A2_FS, k, st_p );
    
    % Possibly do FLD rescaling
    if ( st_p.bp_fld )
        [ a_A2, a_A1, a_A0, alpha ] = fn_ww__util__fld_scaling( a_A2, a_A1, a_A0 );              
    end  
    
    % Do the quadratic eigenvalue calculation
    [ a_w, v_eig ] = polyeig( a_A0, a_A1, a_A2 );
    
    % May have to put the eigenvalues back into the correct scale
    if ( st_p.bp_fld )
        v_eig = v_eig * alpha;
    end
    
    % Now remove eigenvectors corresponding to dunted eigenvalues    
    [ v_eig, a_w ] = fn_ww__calc_re__filter_valid_eigenvalues( v_eig, a_w, st_p );
        
    % Error check
    assert( numel( v_eig ) > 0, 'No valid eigenvalues found after filtering.' );
    
    % Check for critical layer or zero eigenvalue, and process accordingly
    [ max_eig, max_idx ] = max( v_eig );
    if ( max_eig <= U_max || abs( max_eig ) < st_p.fp_zero_tol )
    
        % We're in critical layer terriroty, must be smart
        [ v_c_p(lp_theta), a_c_w_p(:,lp_theta) ] = fn_ww__calc_re__find_valid_eigvecs( st_Dn.v_zm, v_eig, a_w, U_min, U_max );
        
    else 
   
        % We're dominant, so can just use the maximum
        v_c_p(lp_theta) = max_eig;
        
        % Pull eigenvector and normalise (and make sure sign consistent)
        if( nargout > 1 )            
            a_c_w_p( :, lp_theta ) = fn_ww__util__normalise_eigenvectors( a_w( :, max_idx ), 2 );
        end
        
    end
        
    % Undo current shift, if required    
    v_c_p(lp_theta) = v_c_p(lp_theta) - st_p.fp_current_shift;  
    
end
if ( st_p.bp_disp_update ), fprintf( '\n' ); end




end
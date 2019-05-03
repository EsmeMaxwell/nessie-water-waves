function [ v_k_p, a_k_w_p, v_residual ] = fn_ww__calc_re__cl__red_k( st_Dn, v_c, st_r_shear, st_p )
%fn_ww__calc_re__cl__red_k: Calc CL-c disperion relation for reduced profile (solve Rayleigh eqn with free-surface)
% 
%   [ v_k_p, a_k_w_p, v_residual ] = fn_ww__calc_re__cl__red_k( st_Dn, v_c, st_v_shear, st_p )
% 
% Calculates the dispersion relation and eigenvectors. Uses c as parameter
% and k as eigenvalue.
% 
% INPUT
%   st_Dn : Differentiation matrices
%   v_c : Vector of phase velocities
%   st_v_shear : Struct of vectorised shear profile
%   st_p : Parameters
% 
% OUTPUT
%   v_k_p : Vector of k results
%   a_k_w_p : Corresponding eigenvectors, arranged columnwise in k
%   (optional)
%   v_residual : Vector of residuals (optional)
% 
% 
% See also
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__diffmtrx__WR_poldif(),
%   fn_ww__setup__lin_map_Dn_to_mapped(),
%   fn_ww__setup__diffmtrx_mp__WR_chebdif(),
%   fn_ww__setup__shear_fn_to_vec()


% Basic sanity checks
assert( isstruct( st_Dn ) );
assert( isnumeric( v_c ) && 1 == size( v_c, 1 ) );
assert( isstruct( st_r_shear ) );
assert( isstruct( st_p ) );

% Set required mp digits (remember, it doesn't work well in parfor loops)
fn_ww__setup__prepare_mp( st_Dn, st_p );


N = numel( st_Dn.v_z0 );
Nc = numel( v_c );

% U matrices and check crit layer position
if ( st_p.bp_phy_calc )
    
    % Get shear data
    [ st_v_shear, U_nd_min, U_nd_max, U_phy_min, U_phy_max ] = fn_ww__util__get_r_shear_data( st_Dn, st_r_shear, st_p );    
    
    a_U = diag( st_v_shear.v_phy_U + st_p.fp_current_shift );
    a_dU = diag( st_v_shear.v_phy_dU );
    a_ddU = diag( st_v_shear.v_phy_ddU );
    
    U_min = U_phy_min + st_p.fp_current_shift;
    U_max = U_phy_max + st_p.fp_current_shift;    
    
%     a_U = diag( st_v_shear.v_phy_U + st_p.fp_current_shift );
%     a_dU = diag( st_v_shear.v_phy_dU );
%     a_ddU = diag( st_v_shear.v_phy_ddU );
    
%    [ U_min, U_max ] = fn_ww__util__find_U_min_max( st_v_shear.v_phy_U + st_p.fp_current_shift );
else
    
    % Get shear data
    [ st_v_shear, U_nd_min, U_nd_max ] = fn_ww__util__get_r_shear_data( st_Dn, st_r_shear, st_p );
    
    a_U = diag( st_v_shear.v_U + st_p.fp_current_shift );
    a_dU = diag( st_v_shear.v_dU );
    a_ddU = diag( st_v_shear.v_ddU );
    
    U_min = U_nd_min + st_p.fp_current_shift;
    U_max = U_nd_max + st_p.fp_current_shift;
    
%     a_U = diag( st_v_shear.v_U + st_p.fp_current_shift );
%     a_dU = diag( st_v_shear.v_dU );
%     a_ddU = diag( st_v_shear.v_ddU );
    
%    [ U_min, U_max ] = fn_ww__util__find_U_min_max( st_v_shear.v_U + st_p.fp_current_shift );
end
if ( st_p.bp_mp )
    a_I = eye( size( a_U ), 'mp' );
else
    a_I = eye( size( a_U ) );
end




% Precalc
a_B = a_I;
a_B( 1,: ) = 0 * a_B( 1,: );
a_B( end,: ) = 0 * a_B( end,: );



% Allocate memory
if ( st_p.bp_mp )
    v_k_p = zeros( 1, Nc, 'mp' );
else
    v_k_p = zeros( 1, Nc );
end

if ( nargout > 1 )
    if ( st_p.bp_mp )
        a_k_w_p = zeros( N, Nc, 'mp' );        
    else
        a_k_w_p = zeros( N, Nc );
    end
end

if ( nargout > 2 )
    if ( st_p.bp_mp )
        v_residual = zeros( 1, Nc, 'mp' );
    else
        v_residual = zeros( 1, Nc );
    end
end






for lp_c=1:Nc

    if ( st_p.bp_disp_update ), fprintf( '%d ', lp_c ); end
    if ( st_p.bp_disp_update && 0 == mod( lp_c, 20 ) ), fprintf( '\n' ); end
    
    c = v_c(lp_c);
    
    % c dept stuff
    v_q = st_v_shear.v_ddU ./ ( st_v_shear.v_U - c );
    a_Q = diag( v_q );

    % Operator
    a_A = st_Dn.a_D2m - a_Q ;
    a_FS = ( ( a_U - c*a_I )^2 * st_Dn.a_Dm ) - ( ( a_U - c*a_I ) * a_dU + ( 1 / st_p.Fr2 ) * a_I );
    a_A(1,:) = a_FS(1,:);
    a_A(end,: ) = 0 * a_A( end,: );
   
    % Do the eigenvalue calculation
    if ( st_p.bp_ex_btm_bc )

        % Calculate on rows without bottom BC
        [ a_w, a_eig ] = eig( a_A(1:end-1,1:end-1), a_B(1:end-1,1:end-1) );
        
        % Extend back to full length
        a_w = [ a_w; zeros( 1, size( a_w, 2 ) ) ];        
        
    else
        
        % Calc eigenvalue and eigenvector
        [ a_w, a_eig ] = eig( a_A, a_B );
        
    end
        
    v_eigs = diag( a_eig );
    v_eigs( ~isfinite(v_eigs) ) = 0;
    v_eigs( v_eigs < 0 ) = 0;

    % Error check
    assert( numel( v_eigs ) > 0, 'No valid eigenvalues found after filtering.' );
        
    % Get dominant eigenvalue and pull associated eigenvector
    [ mu, max_idx ] = max( v_eigs );
    v_k_p(lp_c) = sqrt( mu );
    
    % Pull eigenvector and normalise (and make sure sign consistent)
    if( nargout > 1 )
        a_k_w_p( :, lp_c ) = fn_ww__util__normalise_eigenvectors( a_w( :, max_idx ), st_p.ip_std_evec_norm  );

        % If we need to calculate residual
        if ( nargout > 2 )            
            v_residual( lp_c ) = norm( ( a_A - mu * a_B ) * a_k_w_p( :, lp_c ), 2 );
        end
    end
    
end
if ( st_p.bp_disp_update ), fprintf( '\n' ); end


end
function [ v_k_CL, a_k_w_CL ] = fn_ww__calc_re__cl__gen_sl_pol_rad_k( st_Dn, v_c, theta, st_v_shear_x, st_v_shear_y, st_p )
%fn_ww__calc_re__cl__gen_sl_pol_rad_k: INCOMPLETE Calc k dispersion relation general radial


error('May need to debug this.')

% Basic sanity checks
assert( isstruct( st_Dn ) );
assert( isnumeric( v_c ) && 1 == size( v_c, 1 ) );
assert( isnumeric( theta ) && 1 = numel( theta ) );
assert( isstruct( st_v_shear_x ) );
assert( isstruct( st_v_shear_y ) );
assert( isstruct( st_p ) );

% Set required mp digits (remember, it doesn't work well in parfor loops)
fn_ww__setup__prepare_mp( st_Dn, st_p );



N = numel( st_Dn.v_z0 );



a_I = eye( size( st_Dn.a_Dm ) );


v_k_CL = zeros( 1, numel(v_c) );
for lp_c=1:numel(v_c)

    c = v_c(lp_c);
    c_x = cos(theta)*c;
    c_y = sin(theta)*c;
    
    % U matrices
    v_U = (1/c) * ( c_x * st_v_shear_x.v_U + c_y * st_v_shear_y.v_U );
    v_dU = (1/c) * ( c_x * st_v_shear_x.v_dU + c_y * st_v_shear_y.v_dU );
    v_ddU = (1/c) * ( c_x * st_v_shear_x.v_ddU + c_y * st_v_shear_y.v_ddU );
    
    % U matrices
    a_U = diag( v_U );
    a_dU = diag( v_dU );
    a_ddU = diag( v_ddU );        

    % Precalc
    a_B = a_I;
    a_B( 1,: ) = 0 * a_B( 1,: );
    a_B( end,: ) = 0 * a_B( end,: );    
    
    % c dept stuff
    v_q = v_ddU ./ ( v_U - c );
    a_Q = diag( v_q );

    % Operator
    a_A = st_Dn.a_D2m - a_Q ;
    a_FS = ( ( a_U - c*a_I )^2 * st_Dn.a_Dm ) - ( ( a_U - c*a_I ) * a_dU + ( 1 / st_p.Fr2 ) * a_I );
    a_A(1,:) = a_FS(1,:);
    a_A(end,: ) = 0 * a_A( end,: );
   
    % Calc eigenvalue and eigen vector
    [ a_w, a_eig ] = eig( a_A, a_B );
    v_eigs = diag( a_eig );
    v_eigs( ~isfinite(v_eigs) ) = 0;
    v_eigs( v_eigs < 0 ) = 0;

    % Error check
    assert( numel( v_eigs ) > 0, 'No valid eigenvalues found after filtering.' );
        
    % Get dominant eigenvalue and pull associated eigenvector
    [ mu, max_idx ] = max( v_eigs );
    v_k_CL(lp_c) = sqrt( mu );
    
    % Pull eigenvector and normalise (and make sure sign consistent)
    if( nargout > 1 )    
        if ( a_w( 1, max_idx ) < 0 )
            a_w( :, max_idx ) = a_w( :, max_idx ) * -1;
        end    
        a_k_w_CL( :, lp_c ) = a_w( :, max_idx ) / norm( a_w( :, max_idx ), 2 );
    end    
    
end


end
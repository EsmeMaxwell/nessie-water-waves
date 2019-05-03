function [ a_k_CL ] = fn_ww__calc_re__cl__red_k__fullspec( st_Dn, v_c, st_v_shear, st_p )
%fn_ww__calc_re__cl__red_k__fullspec: INCOMPLETE Calc dispersion relation for k + full spectrum




% Basic sanity checks
assert( isstruct( st_Dn ) );
assert( isnumeric( v_c ) && 1 == size( v_c, 1 ) );
assert( isstruct( st_v_shear ) );
assert( isstruct( st_p ) );

% Set required mp digits (remember, it doesn't work well in parfor loops)
fn_ww__setup__prepare_mp( st_Dn, st_p );



N = numel( st_Dn.v_z0 );
Nreturn = ceil(N/2);
Nc = numel( v_c );

% U matrices
a_U = diag( st_v_shear.v_U );
a_dU = diag( st_v_shear.v_dU );
a_ddU = diag( st_v_shear.v_ddU );
a_I = eye( size( a_U ) );

% Precalc
a_B = a_I;
a_B( 1,: ) = 0 * a_B( 1,: );
a_B( end,: ) = 0 * a_B( end,: );

a_k_CL = zeros( Nreturn, Nc );

for lp_c=1:Nc

    c = v_c(lp_c);
    
    % c dept stuff
    v_q = st_v_shear.v_ddU ./ ( st_v_shear.v_U - c );
    a_Q = diag( v_q );

    % Operator
    a_A = st_Dn.a_D2m - a_Q ;
    a_FS = ( ( a_U - c*a_I )^2 * st_Dn.a_Dm ) - ( ( a_U - c*a_I ) * a_dU + ( 1 / st_p.Fr2 ) * a_I );
    a_A(1,:) = a_FS(1,:);
    a_A(end,: ) = 0 * a_A( end,: );       
    
    % Calc eigenvalue and eigenvector
    [ v_eigs ] = eig( a_A, a_B );
    
    v_eigs( ~isfinite(v_eigs) ) = [];
    
    v_eigs = sort( v_eigs, 'descend' );
%    v_eigs( v_eigs < 0 ) = 0;

    a_k_CL( :, lp_c ) = v_eigs( 1:Nreturn );
        
end


end
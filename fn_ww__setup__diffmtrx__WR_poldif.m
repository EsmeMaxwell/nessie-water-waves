function [ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( N, b_nst )
%fn_ww__setup__diffmtrx__WR_poldif: Generate diff matrices using WR poldif
%
%   [ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( N, b_nst )
% 
% Constructs the first and second order differentiation matrix, of size
% (N+1), using poldif from the Weideman-Reddy DM suite. The negative sum
% trick can be applied by setting the b_nst flag nonzero. Also returns the
% corresponding z vector of second-kind Chebyshev collocation points on
% [-1,1].
% 
% INPUT
%
% N : Size, actually will return a matrix and vector of size (N+1)
% 
% b_nst : Flag for whether to apply negative sum trick. Anything nonzero
% will enable.
% 
% OUTPUT
%
% st_Dn : struct containing matrices a_D, a_D2, and v_z0.
%
% See also
%   fn_ww__setup__diffmtrx__WR_chebdif(),
%   fn_ww__setup__diffmtrx_mp__WR_chebdif()
%

%v_z0 = cos( pi*(0:N)/N)';

% Ensure float symmetry in the collocation points
v_z0 = sin( pi * ( N - 2*(0:N) ) / ( 2*N) )';

% Call original function asking for diff matrices up to 2nd order
a_DM = fn_ww__ext__diffmtrx__WR_poldif( v_z0, 2 );
a_D = a_DM(:,:,1);
a_D2 = a_DM(:,:,2);



% Negative sum trick (recalculate the diagonal elements to ensure each row
% sum = 0).  For numerical stability, this requires summing the smallest
% elements in each row first

% v_test = ones( N+1, 1 );
% test_orig_D = norm( a_D_nst * v_test )
% test_orig_D2 = norm( a_D2_nst * v_test )
% 
% v_old_D_diag = diag( a_D_nst );
% v_old_D2_diag = diag( a_D2_nst );

if ( b_nst ~= 0 )

    for lp_row=1:(N+1)

        %fprintf( '\n\n%d\n', lp_row );
        D_row_sum = 0;
        D2_row_sum = 0;

        % Create copy of row and eliminate the diagonal element
        v_D_row_no_diag = a_D( lp_row, : );
        v_D_row_no_diag( lp_row ) = 0;
        v_D2_row_no_diag = a_D2( lp_row, : );
        v_D2_row_no_diag( lp_row ) = 0;

        % Get a sort index
        [ ~, v_D_idx ] = sort( abs(v_D_row_no_diag) );
        [ ~, v_D2_idx ] = sort( abs(v_D2_row_no_diag) );

        % Sum, carefully
        for lp_col=1:(N+1)        
            D_row_sum = D_row_sum + v_D_row_no_diag( v_D_idx( lp_col ) );
            D2_row_sum = D2_row_sum + v_D2_row_no_diag( v_D2_idx( lp_col ) );        
        end

        % Replace diagonal entry
        a_D( lp_row, lp_row ) = -D_row_sum;
        a_D2( lp_row, lp_row ) = -D2_row_sum;    

    end

end
    
st_Dn = struct;
st_Dn.v_z0 = v_z0;
st_Dn.a_D = a_D;
st_Dn.a_D2 = a_D2;
st_Dn.N = N;

% 
% 
% diff_D = norm( v_old_D_diag - diag( a_D_nst ) )
% diff_D2 = norm( v_old_D2_diag - diag( a_D2_nst ) )
% 
% test_nst_D = norm( a_D_nst * v_test )
% test_nst_D2 = norm( a_D2_nst * v_test )
% 
% improvement_D = test_orig_D - test_nst_D
% improvement_D2 = test_orig_D2 - test_nst_D2




end
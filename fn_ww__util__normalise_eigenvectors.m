function [ a_w, v_scale ] = fn_ww__util__normalise_eigenvectors( a_w, i_method )
%fn_ww__util__normalise_eigenvectors: Util normalise eigenvectors (columnwise)
%
%   [ a_w, v_scale ] = fn_ww__util__normalise_eigenvectors( a_w, i_method )
%
% Normalises the supplied eigenvectors either by inf, 2-norm, or surface
%
% INPUT
%   a_w : Eigenvectors, arranged columnwise
%   i_method : Method for normalisation
% 
% OUTPUT
%   a_w : Normalised eigenvectors
%   v_scale : Vector of scale factors (when appropriate)
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__calc_re__cl__red_c()


for lp_j=1:size( a_w, 2 )
    
    % Correct by sign so that the surface is always +ve
    if ( a_w( 1, lp_j ) < 0 )
        a_w( :, lp_j ) = a_w( :, lp_j ) * -1;
    end
    
    switch i_method
        case 0  % Normalise against maximum, i.e. inf norm
            v_scale(lp_j) = norm( a_w( :, lp_j ), 'inf' );
            a_w( :, lp_j ) = a_w( :, lp_j ) / v_scale(lp_j);
            max_j = max( a_w( :, lp_j ) );
            min_j = min( a_w( :, lp_j ) );
            if ( abs( max_j ) < abs( min_j ) )
                a_w( :, lp_j ) = a_w( :, lp_j ) * -1;
            end
            continue;
        case 2  % Normalise against 2 norm
            v_scale(lp_j) = norm( a_w( :, lp_j ), 2 );
            a_w( :, lp_j ) = a_w( :, lp_j ) / v_scale(lp_j);
            continue;
        case 100
            assert( abs( a_w( 1, lp_j ) ) > 1e-15, 'Surface near zero, ill-conditioned' );
            v_scale(lp_j) = a_w( 1, lp_j );
            a_w( :, lp_j ) = a_w( :, lp_j ) / v_scale(lp_j);
            continue;
        otherwise
            error( 'No valid method specified' );
    end
        
end


end
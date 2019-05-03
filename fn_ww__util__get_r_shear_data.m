function [ st_v_shear, U_nd_min, U_nd_max, U_phy_min, U_phy_max ] = fn_ww__util__get_r_shear_data( st_Dn, st_r_shear, st_p )
%fn_ww__util__get_r_shear_uminmax: Util get the uminmax data, either from cache or force calculation
% 
%   [ st_crit ] = fn_ww__util__find_U_min_max__fn( st_fn_shear, st_p )
% 
% 
% 

% DEV check, avoid another 2 hours wasted on a subtle flaw...
if ( isfield( st_r_shear, 'v_U' ) || isfield( st_r_shear, 'fn_U' ) )
    error( 'Do not pass st_v_shear or st_fn_shear directly!' );
end


% Do shear discretisation, if necessary (we may already have st_v_shear but
% no st_fn_shear)
if ( isfield( st_r_shear, 'st_fn_shear' ) && ~isfield( st_r_shear, 'st_v_shear' ) )
    [ st_v_shear ] = fn_ww__setup__shear_fn_to_vec( st_Dn, st_r_shear.st_fn_shear, st_p );
else
    % In this case, we should check that the length of existing st_v_shear
    % is corect
    assert( numel( st_Dn.v_z0 ) == numel( st_r_shear.st_v_shear.v_U ), 'Mismatch between st_Dn size and st_v_shear size (or you may need to initialise st_r_shear with fn rather than vec)' );
    st_v_shear = st_r_shear.st_v_shear;
end


% Check if we've already precalculated this
if ( isstruct( st_r_shear.st_crit ) )
    
    % Use precalculated values
    st_crit = st_r_shear.st_crit;

else
   
   % Must calculate, decide whether to do it fn-wise or vec-wise
   if ( st_p.bp_Uminmax_vec || isfield( st_r_shear, 'st_v_shear' ) )              
       fprintf( 'Calculating Uminmax using vector coords\n' );
       [ st_crit ] = fn_ww__util__find_U_min_max__vec( st_Dn, st_v_shear, st_p );
   else
       fprintf( 'Calculating Uminmax using fn optimisation\n' );
       [ st_crit ] = fn_ww__util__find_U_min_max__fn( st_fn_shear, st_p );                    
   end
    
end


U_nd_min = st_crit.U_min;
U_nd_max = st_crit.U_max;

if ( isfield( st_crit, 'U_phy_min' ) )
    U_phy_min = st_crit.U_phy_min;
    U_phy_max = st_crit.U_phy_max;
else
    
    if ( nargout > 3 )
        warning( 'Physical Uminmax not calculated, do not ask for them!' );
    end
    
    U_phy_min = 0;
    U_phy_max = 0;
    
end

end
function [ v_theta_q_rs, v_k_q_rs, v_kx_q_rs, v_ky_q_rs ] = fn_ww__util__create_grid_ordered_pair_polar( v_k_q, v_theta_q, bp_sanitise )
%fn_ww__util__create_grid_ordered_pair_polar: Util forms grid from v_k, v_theta and returns vector of ordered pairs
%
%   [ v_theta_q_rs, v_k_q_rs, v_kx_q_rs, v_ky_q_rs ] = fn_ww__util__create_grid_ordered_pair_polar( v_k_q, v_theta_q, bp_sanitise )
% 
% Forms a grid from the k and theta vector then reshapes and returns a
% vector of ordered pairs.
 

% If we should sanitise the theta vector, i.e. make sure it's on interval
% [-pi,pi)
if ( bp_sanitise )
    [ v_theta_q ] = fn_ww__util__sanitise_theta_query_vec( v_theta_q );
end

% Create grids and then reshape into a vectors
[ a_theta_q, a_k_q ] = meshgrid( v_theta_q, v_k_q );
v_theta_q_rs = reshape( a_theta_q, 1, [] );
v_k_q_rs = reshape( a_k_q, 1, [] );

% If the user wants the corresponding Cartesian coordinates
if ( nargout > 2 )
    [ v_kx_q_rs, v_ky_q_rs ] = pol2cart( v_theta_q_rs, v_k_q_rs );
end




end
function [ v_y, v_dy, v_ddy ] = fn_ww__dev__de_casteljau( v_s, v_q, a, b )
%fn_ww__dev__de_casteljau: DEV Basic implementation of De Casteljau's algorithm
%
%   [ v_y, v_dy, v_ddy ] = fn_ww__dev__de_casteljau( v_s, v_q )
%
% Calculate 1d Bezier curve.  Assumes uniformally spaced v_q on [0,1]. Also
% assumes v_x on [0,1].
%
%



% TODO sort vector input dimensinos, this is silly

if ( size( v_q, 2 ) > 1 )
    v_q = v_q';
end

if ( size( v_s, 1 ) > 1 )
    v_s = v_s';
end


% Remap query pt input vector (warning: we don't need to mess with v_s as
% that's values assumed on uniformly partitioned interval on [0,1])
%v_s = ( v_s - a ) / ( b - a );
v_q = ( v_q - a ) / ( b - a );

N = size( v_s, 2 ) - 1;
Nx = numel( v_q );

a_s = repmat( v_s, Nx, 1 );


for lp_j=1:N
    
    a_n = zeros( size( a_s ) );            
    for lp_i=0:(N-lp_j)
   
%        fprintf( 'j=%d, i=%d\n', lp_j, lp_i );        
        %v_n(lp_i+1) = ( 1 - v_x ) * v_q(lp_i+1) + v_x * v_q(lp_i+2);        
        a_n(:,lp_i+1) = ( 1 - v_q ) .* a_s(:,lp_i+1) + v_q .* a_s(:,lp_i+2);
        
    end

    if ( lp_j > 1 )
        a_dds = a_ds;
    end
    a_ds = a_s;
    a_s = a_n;
    
end

v_y = a_s(:,1);
v_dy = N * ( a_ds(:,2) - a_ds(:,1) );
v_ddy = N * (N-1) * ( a_dds(:,3) - 2 * a_dds(:,2) + a_dds(:,1) );





% General M-dimensional version.  Can't figure out how to vectorise, so
% just limiting to 1d uniform above.

% for lp_j=1:N
%     
%     v_n = zeros( size( v_q ) );            
%     for lp_i=0:(N-lp_j)
%    
%         fprintf( 'j=%d, i=%d\n', lp_j, lp_i );        
%         v_n(:,lp_i+1) = ( 1 - x ) * v_q(:,lp_i+1) + x * v_q(:,lp_i+2);
%         
%     end
% 
%     if ( lp_j > 1 )
%         v_ddq = v_dq;
%     end
%     v_dq = v_q;
%     v_q = v_n;
%     
% end
% 
% y = v_q(:,1);
% dy = N * ( v_dq(:,2) - v_dq(:,1) );
% ddy = N * (N-1) * ( v_ddq(:,3) - 2 * v_ddq(:,2) + v_ddq(:,1) );
% 


end
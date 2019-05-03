function [ a_interp ] = fn_ww__interp__bary_ch_2nd_kind( st_Dn, a_w, st_p, v_q )
%fn_ww__interp__extend_short_eigenvector: 
% 
%   [ v_interp ] = fn_ww__interp__bary_ch_2nd_kind( st_Dn, v_w, st_p, v_q )
% 
% This basically copied from Berrut and Trefethen,
% https://people.maths.ox.ac.uk/trefethen/barycentric.pdf except with
% slight error in the code on p510 fixed.
% 
% TODO : proper credit note
% TODO : optimise
% Note: for now we calculate in true Chebyshev interval [-1,1].
%



% Function 2nd kind evaluation pts, weights.
N = st_Dn.N;
v_x = st_Dn.v_z0;
v_c = [ 1/2; ones( N-1, 1 ); 1/2 ] .* (-1) .^ ((0:N)');

% Extend v_x, v_c to correct columnwise length arrays
Ncols = size( a_w, 2 );
a_f = a_w;
a_x = repmat( v_x, 1, Ncols );
a_c = repmat( v_c, 1, Ncols );

% Map v_q to [-1,1] then extend to array
v_xx = fn_ww__setup__lin_conv_zm_to_z0( v_q, st_p.a, st_p.b );
a_xx = repmat( v_xx, 1, Ncols );

% Do interpolation
a_numer = zeros( size( a_xx ) );
a_denom = zeros( size( a_xx ) );
a_exact = zeros( size( a_xx ) );

for lp_q=1:N+1
    a_xdiff = a_xx - a_x( lp_q, : );
    a_temp = a_c( lp_q, : ) ./ a_xdiff;
    a_numer = a_numer + a_temp .* a_f( lp_q, : );
    a_denom = a_denom + a_temp;
    a_exact( 0 == a_xdiff ) = lp_q;
end

a_interp = a_numer ./ a_denom;

% For now, just loop across the columns... can work out the indexing
% properly later.
for lp_col=1:Ncols
    v_idx = find( a_exact( :, lp_col ) );
    a_interp( v_idx, lp_col ) = a_f( a_exact( v_idx, lp_col ), lp_col );
end




% % Function 2nd kind evaluation pts, weights.
% N = st_Dn.N;
% v_x = st_Dn.v_z0;
% v_c = [ 1/2; ones( N-1, 1 ); 1/2 ] .* (-1) .^ ((0:N)');
% v_f = a_w;
% 
% % Map v_q to [-1,1]
% v_xx = fn_ww__setup__lin_conv_zm_to_z0( v_q, st_p.a, st_p.b );
% 
% % Do interpolation
% v_numer = zeros( size( v_xx ) );
% v_denom = zeros( size( v_xx ) );
% v_exact = zeros( size( v_xx ) );
% 
% for lp_q=1:N+1
%     v_xdiff = v_xx - v_x( lp_q );
%     v_temp = v_c( lp_q ) ./ v_xdiff;
%     v_numer = v_numer + v_temp * v_f( lp_q );
%     v_denom = v_denom + v_temp;
%     v_exact( v_xdiff == 0 ) = lp_q;
% end
% 
% v_interp = v_numer ./ v_denom;
% v_idx = find( v_exact );
% v_interp( v_idx ) = v_f( v_exact( v_idx ) );

end
function [ a_w_dst ] = fn_ww__interp__extend_short_eigenvector( st_Dn_src, a_w_src, st_p_src, v_zm_dst, st_p_dst )
%fn_ww__interp__extend_short_eigenvector: 
% 
%   [ v_w_dst ] = fn_ww__interp__extend_short_eigenvector( st_Dn_src, v_w_src, st_p_src, st_Dn_dst, st_p_dst )
% 
% This extends a function from smaller interval to larger
%
% Note: can do multiple eigenvectors at once but with single domain mapping
%
% Note: assumes ordering
%
% TODO: optimise

% Example: [ v_w_dst ] = fn_ww__interp__extend_short_eigenvector( st_Dn_test, a_c_w_p_test(:,end), st_p_test, st_Dn_160.v_zm, st_p )

% Shorthand
a_src = st_p_src.a;
a_dst = st_p_dst.a;

Ncols = size( a_w_src, 2 );
Ndstrows = size( v_zm_dst, 1 );


% Split dst into the interval we need to inteprolate and the interval which
% should be zero
v_interp_mask = v_zm_dst >= a_src;
v_interp_idx = find( v_interp_mask );

% Ok, calc the required inteprolation point vector
v_q = v_zm_dst( v_interp_idx );
Nqrows = size( v_q, 1 );

% Do the interpolation
[ a_interp ] = fn_ww__interp__bary_ch_2nd_kind( st_Dn_src, a_w_src, st_p_src, v_q );

% Fixup
a_w_dst = zeros( Ndstrows, Ncols );
a_w_dst( 1:Nqrows, : ) = a_interp;




% % Shorthand
% a_src = st_p_src.a;
% a_dst = st_p_dst.a;
% 
% 
% % Split dst into the interval we need to inteprolate and the interval which
% % should be zero
% v_interp_mask = v_zm_dst > a_src;
% v_interp_idx = find( v_interp_mask );
% 
% % Ok, calc the required inteprolation point vector
% v_q = v_zm_dst( v_interp_idx );
% 
% % Do the interpolation
% [ v_interp ] = fn_ww__interp__bary_ch_2nd_kind( st_Dn_src, v_w_src, st_p_src, v_q );
% 
% % Fixup
% v_w_dst = zeros( size( v_zm_dst ) );
% v_w_dst(1:numel(v_interp)) = v_interp;

end
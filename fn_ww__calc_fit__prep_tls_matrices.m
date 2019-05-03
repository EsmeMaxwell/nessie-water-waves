function [ st_fit_tlscb_mtrx ] = fn_ww__calc_fit__prep_tls_matrices( st_fit_cb, tls_mode, tls_alpha )
%fn_ww__calc_fit__prep_tls_matrices: Calc TLS Chebyshev precalc matrices
%
%   [ st_fit_tlscb_mtrx ] = fn_ww__calc_fit__prep_tls_matrices( st_fit_cb, tls_mode, tls_alpha )
%
% From CB basis, calculates the necessary Tikhonov SVD matrices for
% regularised least squares calc
%
% TAGS: WWERRINSHEAR
% 
% See also
%   fn_ww__calc_fit__prep_cb_basis(),
%   fn_ww__calc_fit__ols(),
%   fn_ww__calc_fit__prep_lin_surf_extrapolate()


st_fit_tlscb_mtrx = struct;

Mcb = size( st_fit_cb.a_T, 2 );

switch tls_mode
    case 0
        a_L = eye( Mcb );
        rows_T = size( a_L, 1 );
    case 2
        a_L = zeros( Mcb-2, Mcb );
        rows_T = size( a_L, 1 );
        for lp_i=1:rows_T
            a_L(lp_i,lp_i:(lp_i+2)) = [ 1 -1 2 ];
        end
    otherwise
        error( 'Error' );     
end


% Tikhonov using SVD
a_Tk = [ st_fit_cb.a_T; tls_alpha * a_L ];
[ a_U, a_S, a_V ] = svd( a_Tk );
a_Spinv = ( a_S )';
a_Spinv( a_Spinv ~= 0 ) = 1 ./ a_Spinv( a_Spinv ~= 0 );  % Eugh! Fix me.


st_fit_tlscb_mtrx.tls_mode = tls_mode;
st_fit_tlscb_mtrx.tls_alpha = tls_alpha;
st_fit_tlscb_mtrx.rows_T = rows_T;
st_fit_tlscb_mtrx.a_L = a_L;
st_fit_tlscb_mtrx.a_U = a_U;
st_fit_tlscb_mtrx.a_S = a_S;
st_fit_tlscb_mtrx.a_V = a_V;
st_fit_tlscb_mtrx.a_Spinv = a_Spinv;




end
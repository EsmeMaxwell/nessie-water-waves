function [ v_k ] = fn_ww__util__create_k_vec( k_min, k_max, Nk, ip_mode, fp_safety_margin )
%fn_ww__util__create_k_vec: Util creates k vectors
%
%   [ v_k ] = fn_ww__util__create_k_vec( k_min, k_max, Nk, ip_mode, fp_safety_margin )
%
% Generates k vectors in a linear, logarithmic, mixed, linear-log or
% quadratic arrangement. Can also generate random vectors in the same
% distributions. This function should be used throughout WW_LIB for
% consistency.
%
% INPUT
%   k_min, k_max : k interval to generate vector within
%   Nk : Number of elements to include in k vector
%   ip_mode : Type of distribution to use, integer id as follows,
%       0=linear
%       1=log
%       2=mix of linear and log
%       3=quadtratic
%       10=random (from lin dist)
%       11=random (from log dist)
%       12=random (from mixed dist)
%       13=random (from quad dist)
%   fp_safety_margin : Contract interval at each end by this amount. Should
% be used when generating vectors being supplied to the PF interpolation
% functions.
%
% OUTPUT
%   v_k : the expected k vector
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__calc_re__cl__red_c(),
%   fn_ww__util__create_theta_vec(),
%   fn_ww__calc_re__cl__gen_sl_pol_rad_c()


% Error checks
assert( k_min >= 0, 'k_min should be non-negative.' );
assert( k_max >= 0, 'k_max should be non-negative.' );
assert( k_min < k_max, 'k_min should be strictly less than k_max.' );
assert( Nk > 2, 'Nk should be greater than 2.' );
assert( fp_safety_margin >= 0, 'Safety margin should be non-negative.' );
assert( fp_safety_margin < ( k_max - k_min ), 'Safety margin cannot push values outside of computation range.' );

if ( 2 == ip_mode ), warning( 'Need to fix type 2, have duplicates' ); end

% Create initial safety margin
k_min_s = k_min + fp_safety_margin;
k_max_s = k_max - fp_safety_margin;

fn_poly_space_lin = @( x1, x2, deg, n ) x1 + ( x2 - x1 ) * linspace( 0, 1, n ).^deg;
fn_poly_space_rand = @( x1, x2, deg, n ) x1 + ( x2 - x1 ) * rand( 1, n ).^deg;


switch ip_mode
    case 0
        v_k = linspace( k_min_s, k_max_s, Nk );
    case 1
        v_k = logspace( log10( k_min_s ), log10( k_max_s ), Nk );
    case 2
        v_k_lin = linspace( k_min_s, k_max_s, floor(Nk/2) );
        v_k_log = logspace( log10( k_min_s ), log10( k_max_s ), ceil(Nk/2) );
        v_k = sort( [ v_k_lin v_k_log ] );
    case 3 
        v_k = fn_poly_space_lin( k_min_s, k_max_s, 2, Nk );
    case 10
        v_k = k_min_s + ( k_max_s - k_min_s ) * rand( 1, Nk );
        v_k = sort( v_k );
    case 11
        v_k = log10( k_min_s ) + ( log10( k_max_s ) - log10( k_min_s ) ) * rand( 1, Nk );
        v_k = 10 .^ v_k;
        v_k = sort( v_k );
    case 12
        v_k_lin = k_min_s + ( k_max_s - k_min_s ) * rand( 1, floor(Nk/2) );
        v_k_log = log10( k_min_s ) + ( log10( k_max_s ) - log10( k_min_s ) ) * rand( 1, ceil(Nk/2) );
        v_k_log = 10 .^ v_k_log;
        v_k = sort( [ v_k_lin v_k_log ] );
    case 13
        v_k = fn_poly_space_rand( k_min_s, k_max_s, 2, Nk );
        v_k = sort( v_k );
    otherwise
        error( 'Must select a valid mode for k range generation.' );
end


% Check we're still inside the correct ranges
v_k( v_k < k_min_s ) = k_min_s;
v_k( v_k > k_max_s ) = k_max_s;


end
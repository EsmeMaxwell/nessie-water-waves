function [ st_fit_cb ] = fn_ww__calc_fit__prep_cb_basis( poly_deg, v_zs, hs )
%fn_ww__calc_fit__prep_cb_basis: Calc Chebyshev basis matrices used in OLS and TLS fits
% 
%   [ st_fit_cb ] = fn_ww__calc_fit__prep_cb_basis( poly_deg, v_zs, hs )
%  
% Takes in a vector of depth abscissas, remaps to Chebyshev interval and
% calculates the basis matrix.
%
% The OLS and TLS methods can do extrapolation to a small extent but it's
% dangerous, i.e. you can supply a smaller interval here compared to the
% subsequent interpolation interval. May however be better to use the
% linear extrapolation function and do poly interp on that reconstructed
% interval.
%
% Note if hs is supplied then we're creating a cb basis for interpolation
% and it's going to go past the standard cb domain, i.e. we take -hs as our
% maximum.
%
% INPUT
%
%   poly_deg : Degree of polynomial to use
%
%   v_zs : z sample locations
%
% TAGS: WWERRINSHEAR
%
% See also
%   fn_ww__calc_fit__prep_lin_surf_extrapolate(),
%   fn_ww__calc_fit__prep_ols_matrices(),
%   fn_ww__calc_fit__prep_tls_matrices(),
%   fn_ww__calc_fit__ols(),
%   fn_ww__calc_fit__tls()




% The vector must be in the correct order
assert( isnumeric( v_zs ) && 1 == size( v_zs, 2 ), 'z sample location vector must be column vector' );
assert( v_zs(1) > v_zs(end), 'z sample location vector must be in descending order' );

% Fix the off-by-one problem
Mcb = poly_deg + 1;
Nz = numel( v_zs );
st_fit_cb = struct;

% Remap to [-1,1] ... will fix this properly later
min_v_zm = min(v_zs);
max_v_zm = max(v_zs);
if ( nargin > 2 )
    max_v_zm = -hs;
end
v_z0 = ( ( v_zs - min_v_zm ) - ( max_v_zm - v_zs) ) / ( max_v_zm - min_v_zm );

% Chebyshev 1st kind
a_T = ones( Nz, Mcb );
a_T(:,2) = v_z0;
for lp_k=3:Mcb
    a_T(:,lp_k) = 2 * v_z0 .* a_T(:,lp_k-1) - a_T(:,lp_k-2);
end

% Chebyshev 2nd kind
a_U = ones(Nz,Mcb);
a_U(:,2) = 2*v_z0;
for lp_k=3:Mcb    
    a_U(:,lp_k) = 2 * v_z0 .* a_U(:,lp_k-1) - a_U(:,lp_k-2);
end

% Chebyshev 1st derivative first kind
a_dT = zeros(Nz,Mcb);
a_dT(:,2) = a_U(:,1);
for lp_k=3:Mcb
    a_dT(:,lp_k) = (lp_k-1) * a_U(:,lp_k-1);
end

% Chebyshev 2nd derivative first kind
a_ddT = zeros(Nz,Mcb);
for lp_k=2:Mcb
    a_ddT(:,lp_k) = (lp_k-1) * ( lp_k * a_T(:,lp_k) - a_U(:,lp_k) ) ./ ( v_z0.^2 - 1  );
    
    % These are to avoid divide-by-zero problems.  The top one may not be
    % required if using hs
    if ( abs( v_z0(1) - 1 ) < 1e-14 ) 
        a_ddT(1,lp_k) = ( (lp_k-1)^4 - (lp_k-1)^2 ) / 3;
    end
    a_ddT(end,lp_k) = (-1)^(lp_k-1) * ( (lp_k-1)^4 - (lp_k-1)^2 ) / 3;           
end


st_fit_cb.min_v_zm = min_v_zm;
st_fit_cb.max_v_zm = max_v_zm;
st_fit_cb.v_z0 = v_z0;
st_fit_cb.v_zs = v_zs;
st_fit_cb.a_T = a_T;
st_fit_cb.a_dT = a_dT;
st_fit_cb.a_ddT = a_ddT;
st_fit_cb.a_U = a_U;

end
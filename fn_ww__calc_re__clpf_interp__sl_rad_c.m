function [ a_y0_PF_interp, a_dy0_PF_interp ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp, v_k_q, b_vector )
%fn_ww__calc_re__clpf_interp__sl_rad_c: Calc interpolant for PF-R-r
%
%   [ a_y0_PF_interp, a_dy0_PF_interp ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp, v_k_q, b_vector )
%
% Optimised Hermite interpolation code for path-following algorithm in
% radial direction.
% 
% INPUT
%   st_dp : Struct of control point data
%   v_k_q : Vector of k query points
%   b_vector : Flag, true if eigenvector requested; false if only
%   eigenvalue
% 
% OUTPUT
%   a_y0_PF_interp : Eigenvector, eigenvalues (last row).
%   a_dy0_PF_interp : Derivatives of eigenvector, eigenvalues (last row).
% 
% TAGS: CORE, SISCPFLIB
%
% COPYRIGHT
%   Based on code by S. Loisel / P. Maxwell 2016
%
% See also
%   fn_ww__calc_re__clpf_dp__red_c()



% Note: output is row vector to match with std, which outputs row vector to
% match with eigenvector configuration


q_min = min( [ v_k_q(1) v_k_q(end) ] );
q_max = max( [ v_k_q(1) v_k_q(end) ] );
ctl_min = min( [ st_dp.v_k(1) st_dp.v_k(end) ] );
ctl_max = max( [ st_dp.v_k(1) st_dp.v_k(end) ] );
assert( q_min >= ctl_min && q_max <= ctl_max, 'Interpolation query points not in range of control points!' );


if( st_dp.v_k(1) > st_dp.v_k(end))
    [~,~,v_k_idxs] = histcounts(-v_k_q,-st_dp.v_k);
else
    [~,~,v_k_idxs] = histcounts(v_k_q,st_dp.v_k);
end


if ( 1 == b_vector ) 
    % Shouldn'y have to reassemble but, meh
    a_y = [ st_dp.a_w; st_dp.v_c ];
    a_dy = [ st_dp.a_dw; st_dp.v_dc ];
    a_ymid = [ st_dp.a_wmid; st_dp.v_cmid ];
else
    a_y = st_dp.v_c;
    a_dy = st_dp.v_dc;
    a_ymid = st_dp.v_cmid;
end    

a_y0_PF_interp = zeros(size(a_y,1),length(v_k_idxs));
a_dy0_PF_interp = zeros(size(a_y,1),length(v_k_idxs));

v_k1 = st_dp.v_k( v_k_idxs );
v_k2 = st_dp.v_k( v_k_idxs+1 );
v_km = ( v_k1 + v_k2 ) / 2;
v_dk = v_k2 - v_km;

v_s = ( v_k_q - v_km ) ./ v_dk;   % must be a row vector

a_bf = [-1/4*v_s.*(2*v_s+3).*(v_s-1).^2
          -1/4*v_s.*(1+v_s).*(v_s-1).^2
          (v_s-1).^2.*(v_s+1).^2
          -1/4*v_s.*(2*v_s-3).*(v_s+1).^2
          -1/4*v_s.*(1-v_s).*(v_s+1).^2];
      

a_precalc_y1 = a_y(:,v_k_idxs) .* a_bf(1,:);  % elementwise multiplication by row!
a_precalc_dy1 = a_dy(:,v_k_idxs) .* a_bf(2,:);
a_precalc_ymid = a_ymid(:,v_k_idxs) .* a_bf(3,:);
a_precalc_y2 = a_y(:,v_k_idxs+1) .* a_bf(4,:);
a_precalc_dy2 = a_dy(:,v_k_idxs+1) .* a_bf(5,:);

%a_y0_PF_interp = sum( a_elementwise );
a_y0_PF_interp = a_precalc_y1 + a_precalc_dy1 .* v_dk + a_precalc_ymid + a_precalc_y2 + a_precalc_dy2 .* v_dk;
    
% Process derivative too
if(nargout>1)
	a_dbf = [-(1/4*(v_s-1)).*(8*v_s-3).*(1+v_s)
               -(1/4*(v_s-1)).*(4*v_s.^2+v_s-1)
               4*v_s.*(1+v_s).*(v_s-1)
               -(1/4*(v_s-1)).*(8*v_s+3).*(1+v_s)
               (1/4*(1+v_s)).*(4*v_s.^2-v_s-1)];
           
    a_dy_precalc_y1 = a_y(:,v_k_idxs) .* a_dbf(1,:);  % elementwise multiplication by row!
    a_dy_precalc_dy1 = a_dy(:,v_k_idxs) .* a_dbf(2,:);
    a_dy_precalc_ymid = a_ymid(:,v_k_idxs) .* a_dbf(3,:);
    a_dy_precalc_y2 = a_y(:,v_k_idxs+1) .* a_dbf(4,:);
    a_dy_precalc_dy2 = a_dy(:,v_k_idxs+1) .* a_dbf(5,:);           
           
	a_dy0_PF_interp = a_dy_precalc_y1 ./ v_dk + a_dy_precalc_dy1 + a_dy_precalc_ymid ./ v_dk + a_dy_precalc_y2 ./ v_dk + a_dy_precalc_dy2;
    
end
    
    
return
% 
% 
% for lp_k=1:length(v_k_idxs)
%     k_idx = v_k_idxs(lp_k);
%     k1 = st_dp.v_k(k_idx);
%     k2 = st_dp.v_k(k_idx+1);
%     km = (k1+k2)/2;
%     dk = k2-km;
%     c1 = a_y(:,k_idx);
%     c2 = a_y(:,k_idx+1);
%     dc1 = a_dy(:,k_idx);
%     dc2 = a_dy(:,k_idx+1);
%     cm = a_ymid(:,k_idx);
%     s = (v_k(lp_k)-km)/dk;
%     bf = [-1/4*s.*(2*s+3).*(s-1).^2
%           -1/4*s.*(1+s).*(s-1).^2
%           (s-1).^2.*(s+1).^2
%           -1/4*s.*(2*s-3).*(s+1).^2
%           -1/4*s.*(1-s).*(s+1).^2]
%     a_y0_PF_interp(:,lp_k) = [c1 dk*dc1 cm c2 dk*dc2]*bf;
%     
%     % Derivative
%     if(nargout>1)
%         dbf = [-(1/4*(s-1)).*(8*s-3).*(1+s)
%                -(1/4*(s-1)).*(4*s.^2+s-1)
%                4*s*(1+s).*(s-1)
%                -(1/4*(s-1)).*(8*s+3).*(1+s)
%                (1/4*(1+s)).*(4*s.^2-s-1)];
%         a_dy0_PF_interp(:,lp_k) = [c1./dk dc1 cm./dk c2./dk dc2]*dbf;
%     end    
% end    



end
function [ a_y0_PF_interp, a_dy0_PF_interp ] = fn_ww__calc_re__clpf_interp__sl_ang_c( st_dp_angular, v_theta_q, b_vector )
%fn_ww__calc_re__clpf_interp__sl_ang_c: Calc interpolant for PF-R-a
%
%   [ a_y0_PF_interp, a_dy0_PF_interp ] = fn_ww__calc_re__clpf_interp__sl_ang_c( st_dp_angular, v_theta_q, b_vector )
%
% Optimised Hermite interpolation code for path-following algorithm
% in angular direction.
% 
% INPUT
%   st_dp_angular : Struct of control point data
%   v_theta_q : Vector of theta query points
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


% q_min = min( [ v_theta(1) v_theta(end) ] );
% q_max = max( [ v_theta(1) v_theta(end) ] );
% ctl_min = min( [ st_dp_angular.v_theta(1) st_dp_angular.v_theta(end) ] );
% ctl_max = max( [ st_dp_angular.v_theta(1) st_dp_angular.v_theta(end) ] );
% assert( q_min >= ctl_min && q_max <= ctl_max, 'Interpolation query points not in range of control points!' );

% Check control points are correctly arranged
assert( abs( st_dp_angular.v_theta(1) ) == pi && abs( st_dp_angular.v_theta(end) ) == pi, 'Control points must be on full [-pi,pi] domain.' );

% Sanitise the theta query vector to be in [-pi,pi)
[ v_theta_q ] = fn_ww__util__sanitise_theta_query_vec( v_theta_q );

if( st_dp_angular.v_theta(1) > st_dp_angular.v_theta(end))
    [~,~,v_theta_idxs] = histcounts(-v_theta_q,-st_dp_angular.v_theta);
else
    [~,~,v_theta_idxs] = histcounts(v_theta_q,st_dp_angular.v_theta);
end


if ( 1 == b_vector ) 
    % Shouldn'y have to reassembly but, meh
    a_y = [ st_dp_angular.a_w; st_dp_angular.v_c ];
    a_dy = [ st_dp_angular.a_dw; st_dp_angular.v_dc ];
    a_ymid = [ st_dp_angular.a_wmid; st_dp_angular.v_cmid ];
else
    a_y = st_dp_angular.v_c;
    a_dy = st_dp_angular.v_dc;
    a_ymid = st_dp_angular.v_cmid;
end    

a_y0_PF_interp = zeros(size(a_y,1),length(v_theta_idxs));
a_dy0_PF_interp = zeros(size(a_y,1),length(v_theta_idxs));







v_theta1 = st_dp_angular.v_theta( v_theta_idxs );
v_theta2 = st_dp_angular.v_theta( v_theta_idxs+1 );
v_thetam = ( v_theta1 + v_theta2 ) / 2;
v_dtheta = v_theta2 - v_thetam;

v_s = ( v_theta_q - v_thetam ) ./ v_dtheta;   % must be a row vector

a_bf = [-1/4*v_s.*(2*v_s+3).*(v_s-1).^2
          -1/4*v_s.*(1+v_s).*(v_s-1).^2
          (v_s-1).^2.*(v_s+1).^2
          -1/4*v_s.*(2*v_s-3).*(v_s+1).^2
          -1/4*v_s.*(1-v_s).*(v_s+1).^2];
      

a_precalc_y1 = a_y(:,v_theta_idxs) .* a_bf(1,:);  % elementwise multiplication by row!
a_precalc_dy1 = a_dy(:,v_theta_idxs) .* a_bf(2,:);
a_precalc_ymid = a_ymid(:,v_theta_idxs) .* a_bf(3,:);
a_precalc_y2 = a_y(:,v_theta_idxs+1) .* a_bf(4,:);
a_precalc_dy2 = a_dy(:,v_theta_idxs+1) .* a_bf(5,:);

%a_y0_PF_interp = sum( a_elementwise );
a_y0_PF_interp = a_precalc_y1 + a_precalc_dy1 .* v_dtheta + a_precalc_ymid + a_precalc_y2 + a_precalc_dy2 .* v_dtheta;
    
% Process derivative too
if(nargout>1)
	a_dbf = [-(1/4*(v_s-1)).*(8*v_s-3).*(1+v_s)
               -(1/4*(v_s-1)).*(4*v_s.^2+v_s-1)
               4*v_s.*(1+v_s).*(v_s-1)
               -(1/4*(v_s-1)).*(8*v_s+3).*(1+v_s)
               (1/4*(1+v_s)).*(4*v_s.^2-v_s-1)];
           
    a_dy_precalc_y1 = a_y(:,v_theta_idxs) .* a_dbf(1,:);  % elementwise multiplication by row!
    a_dy_precalc_dy1 = a_dy(:,v_theta_idxs) .* a_dbf(2,:);
    a_dy_precalc_ymid = a_ymid(:,v_theta_idxs) .* a_dbf(3,:);
    a_dy_precalc_y2 = a_y(:,v_theta_idxs+1) .* a_dbf(4,:);
    a_dy_precalc_dy2 = a_dy(:,v_theta_idxs+1) .* a_dbf(5,:);           
           
	a_dy0_PF_interp = a_dy_precalc_y1 ./ v_dtheta + a_dy_precalc_dy1 + a_dy_precalc_ymid ./ v_dtheta + a_dy_precalc_y2 ./ v_dtheta + a_dy_precalc_dy2;
    
end
    
    
return


% 
% 
% 
% for j=1:length(P)
%     p = P(j);
%     theta1 = st_dp_angular.v_theta(p);
%     theta2 = st_dp_angular.v_theta(p+1);
%     thetam = (theta1+theta2)/2;
%     dtheta = theta2-thetam;
%     c1 = a_y(:,p);
%     c2 = a_y(:,p+1);
%     dc1 = a_dy(:,p);
%     dc2 = a_dy(:,p+1);
%     cm = a_ymid(:,p);
%     s = (v_theta(j)-thetam)/dtheta;
%     bf = [-1/4*s.*(2*s+3).*(s-1).^2
%           -1/4*s.*(1+s).*(s-1).^2
%           (s-1).^2.*(s+1).^2
%           -1/4*s.*(2*s-3).*(s+1).^2
%           -1/4*s.*(1-s).*(s+1).^2];
%     a_y0_PF_interp(:,j) = [c1 dtheta*dc1 cm c2 dtheta*dc2]*bf;
%     
%     % Derivative
%     if(nargout>1)
%         dbf = [-(1/4*(s-1)).*(8*s-3).*(1+s)
%                -(1/4*(s-1)).*(4*s.^2+s-1)
%                4*s*(1+s).*(s-1)
%                -(1/4*(s-1)).*(8*s+3).*(1+s)
%                (1/4*(1+s)).*(4*s.^2-s-1)];
%         a_dy0_PF_interp(:,j) = [c1./dtheta dc1 cm./dtheta c2./dtheta dc2]*dbf;
%     end        
% end    


end
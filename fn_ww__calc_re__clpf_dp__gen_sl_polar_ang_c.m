function [ st_dp ] = fn_ww__calc_re__clpf_dp__gen_sl_polar_ang_c( st_Dn, k, st_v_shear_x, st_v_shear_y, tol, st_p )
%fn_ww__calc_re__clpf_dp__gen_sl_polar_ang_c: INCOMPLETE Calc PF-c disperion relation for angular radial profile
%
%   [ st_dp ] = fn_ww__calc_re__clpf_dp__gen_sl_polar_ang_c( st_Dn, k, st_v_shear_x, st_v_shear_y, tol, st_p )
%
% Calculates positive branch of phase velocity for given wave number and
% shear profile using path-following method.
%
% Polar version, calculates in angular curves. This is NOT!! of the same
% ilk as the k->c or c->k versions. This version must use a theta dependent
% shear profile.
%
% Based on code by S. Loisel & P. Maxwell
%
% INPUT
%   st_Dn : Collocation struct, obtain from fn_CL__setup__CL_poldif_nst
%   k_min, k_max : Range of k to calculate within
%   v_zm : Collocation z vector _mapped_, obtained from fn_CL__setup__lin_conv_z0_to_zm
%   st_v_shear : Shear profile sampled at collocation points, use fn_ww__setup__shear_fn_to_vec
%   k0 : Starting point, make it in the middle of the domain, usualy k=2 is good
%   tol : Tolerance for numerical integration, 1e-10 is nice
%
% OUTPUT
%   st_dp : Struct of output that should then be fed into fn_ww__calc_re__clpf_interp__sl_rad_c
%


warning('TODO May want to move the CL solve form here out to a proper function file' );

bp_pf_use_eig_mp = false;

% Basic sanity checks
assert( isstruct( st_Dn ) );
assert( isstruct( st_v_shear_x ) );
assert( isstruct( st_v_shear_y ) );
assert( isstruct( st_p ) );

% Process the control parameters
if ( isfield( st_p, 'bp_pf_use_eig_mp' ) && st_p.bp_pf_use_eig_mp )
    bp_pf_use_eig_mp = st_p.bp_pf_use_eig_mp;    
    st_Dn = st_p.st_Dn_mp;
end




a_D2m = st_Dn.a_D2m;
a_Dm = st_Dn.a_Dm;

N = numel( st_Dn.v_z0 );

a_I = eye( size( a_Dm ) );
a_0I = 0 * a_I;




% Setup numerical integration material (Butcher tableau)
a_BT = [    0    0          0           0          0         0             0        0
            1/5  1/5        0           0          0         0             0        0
            3/10 3/40       9/40        0          0         0             0        0
            4/5  44/45      -56/15      32/9       0         0             0        0
            8/9  19372/6561 -25360/2187 64448/6561 -212/729  0             0        0
            1    9017/3168  -355/33     46732/5247 49/176    -5103/18656   0        0
            1    35/384     0           500/1113   125/192   -2187/6784    11/84    0
            0    35/384     0           500/1113   125/192   -2187/6784    11/84    0
            0    5179/57600 0           7571/16695 393/640   -92097/339200 187/2100 1/40 ];

v_BT_c = a_BT(1:7,1);
a_BT_Z = ( a_BT(1:7,2:end) )';
v_BT_b = ( a_BT(8,2:end) )';
v_BT_b_err = ( a_BT(9,2:end) )' - v_BT_b;
v_BT_b_mid = [  6025192743/30085553152
                0
                51252292925/65400821598
                -2691868925/45128329728
                187940372067/1594534317056
                -1776094331/19743644256
                11237099/235043384];




% Precalc
a_B = a_I;
a_B( 1,: ) = 0 * a_B( 1,: );
a_B( end,: ) = 0 * a_B( end,: );

a_R = a_I;
a_R(1,1) = 0;
a_R(end,end) = 0;
a_R_short = a_R(1:end-1,1:end-1);




% c_max = 1.5;
% c_min = 0.3;
%c0 = 1.0;
[ v_y0 ] = fn__CL_solve( 0 );

% Calculate the DP integral points
[ st_sol_fw ] = fn__dormand_prince( 0, pi, v_y0, tol );
[ st_sol_bw ] = fn__dormand_prince( 0, -pi, v_y0, tol );

% Process the results from the forward and backward integration
v_c_PF_fw_p = fliplr( st_sol_fw.a_y( end,: ) );
v_c_PF_bw_p = st_sol_bw.a_y( end,: );

a_w_PF_fw_p = fliplr( st_sol_fw.a_y( 1:end-1,: ) );
a_w_PF_bw_p = st_sol_bw.a_y( 1:end-1,: );

v_dc_PF_fw_p = fliplr( st_sol_fw.a_dy( end,: ) );
v_dc_PF_bw_p = st_sol_bw.a_dy( end,: );

a_dw_PF_fw_p = fliplr( st_sol_fw.a_dy( 1:end-1,: ) );
a_dw_PF_bw_p = st_sol_bw.a_dy( 1:end-1,: );

v_cmid_PF_fw_p = fliplr( st_sol_fw.a_ymid( end,: ) );
v_cmid_PF_bw_p = st_sol_bw.a_ymid( end,: );

a_wmid_PF_fw_p = fliplr( st_sol_fw.a_ymid( 1:end-1,: ) );
a_wmid_PF_bw_p = st_sol_bw.a_ymid( 1:end-1,: );

v_theta_PF_fw = fliplr( st_sol_fw.v_t );
v_theta_PF_bw = st_sol_bw.v_t;


% Because the last entry of fw == first entry of bw, we need to delete one.
% But first we take a mean, just incase there is a discrepancy.  We do not
% however do this for the mid values, obviously.
v_c_PF_bw_p(1) = 0.5 * ( v_c_PF_fw_p(end) + v_c_PF_bw_p(1) );
a_w_PF_bw_p(:,1) = 0.5 * ( a_w_PF_fw_p(:,end) + a_w_PF_bw_p(:,1) );

v_dc_PF_bw_p(1) = 0.5 * ( v_dc_PF_fw_p(end) + v_dc_PF_bw_p(1) );
a_dw_PF_bw_p(:,1) = 0.5 * ( a_dw_PF_fw_p(:,end) + a_dw_PF_bw_p(:,1) );

v_theta_PF_bw(1) = 0.5 * ( v_theta_PF_fw(end) + v_theta_PF_bw(1) );

v_c_PF_fw_p(end) = [];
a_w_PF_fw_p(:,end) = [];

v_dc_PF_fw_p(end) = [];
a_dw_PF_fw_p(:,end) = [];

v_theta_PF_fw(end) = [];


% Update counter
feval_tot = st_sol_fw.feval + st_sol_bw.feval;
adjeval_tot = st_sol_fw.adjeval + st_sol_bw.adjeval;

% Finally concatenate results
v_c_PF_DP_p = [ v_c_PF_fw_p v_c_PF_bw_p ];
a_w_PF_DP_p = [ a_w_PF_fw_p a_w_PF_bw_p ];

v_dc_PF_DP_p = [ v_dc_PF_fw_p v_dc_PF_bw_p ];
a_dw_PF_DP_p = [ a_dw_PF_fw_p a_dw_PF_bw_p ];

v_cmid_PF_DP_p = [ v_cmid_PF_fw_p v_cmid_PF_bw_p ];
a_wmid_PF_DP_p = [ a_wmid_PF_fw_p a_wmid_PF_bw_p ];

v_theta_PF_DP = [ v_theta_PF_fw v_theta_PF_bw ];


st_dp = struct;
st_dp.k = k;
st_dp.tol = tol;
st_dp.v_theta = v_theta_PF_DP;
st_dp.v_c = v_c_PF_DP_p;
st_dp.a_w = a_w_PF_DP_p;
st_dp.v_dc = v_dc_PF_DP_p;
st_dp.a_dw = a_dw_PF_DP_p;
st_dp.v_cmid = v_cmid_PF_DP_p;
st_dp.a_wmid = a_wmid_PF_DP_p;
st_dp.feval_tot = feval_tot;
st_dp.adjeval_tot = adjeval_tot;







% This seems to work ok...
function [ v_y ] = fn__CL_solve( theta )
   
    % U matrices
    v_U = cos(theta) * st_v_shear_x.v_U + sin(theta) * st_v_shear_y.v_U;
    v_dU = cos(theta) * st_v_shear_x.v_dU + sin(theta) * st_v_shear_y.v_dU;
    v_ddU = cos(theta) * st_v_shear_x.v_ddU + sin(theta) * st_v_shear_y.v_ddU;
    a_U = diag( v_U );
    a_dU = diag( v_dU );
    a_ddU = diag( v_ddU );        

    % Precomputations
    a_A0_precomp = (  a_U * a_D2m - a_ddU  );
    a_A2 = 0 * a_I;

    % FS
    a_A0_FS = ( a_U^2 * a_Dm - a_dU * a_U - ( 1 / st_p.Fr2 ) * a_I );
    a_A1_FS = ( -2 * a_U * a_Dm + a_dU );
    a_A2_FS = a_Dm;
        
    % Main equations
    a_A0 = ( a_A0_precomp - k^2 * a_U );
    a_A1 = -a_D2m + k^2 * a_I;

    a_A0(1,:) = a_A0_FS(1,:);
    a_A1(1,:) = a_A1_FS(1,:);
    a_A2(1,:) = a_A2_FS(1,:);

    a_A0(end,:) = 0 * a_A0(end,:);
    a_A1(end,:) = 0 * a_A1(end,:);
    a_A2(end,:) = 0 * a_A2(end,:);
        
    [ a_w, v_eigs ] = polyeig( a_A0, a_A1, a_A2 );
%    v_eigs = diag( a_eig );
    v_eigs( ~isfinite(v_eigs) ) = 0;
    v_eigs( v_eigs < 0 ) = 0;
    
    % Get max eigenvalue and pull associated eigenvector
    [ c, c_idx ] = max( v_eigs );
    v_w = a_w( :, c_idx );
%    v_y = [ v_w(2:end-1); mu ];   % strip off the first and last
    v_y = [ v_w; c ];   % strip off the first and last
    
end






% This needs to change to a theta dependence
function [ v_dy ] = fn__ode_f( theta, v_y )
   
    %v_w = v_y(1:end-1);
    v_w = v_y(1:end-2);    
    c = v_y(end);
    
    % U matrices
    v_U = cos(theta) * st_v_shear_x.v_U + sin(theta) * st_v_shear_y.v_U;
    v_dU = cos(theta) * st_v_shear_x.v_dU + sin(theta) * st_v_shear_y.v_dU;
    v_ddU = cos(theta) * st_v_shear_x.v_ddU + sin(theta) * st_v_shear_y.v_ddU;
    
    v_dtheta_U = -sin(theta) * st_v_shear_x.v_U + cos(theta) * st_v_shear_y.v_U;
    v_dtheta_dU = -sin(theta) * st_v_shear_x.v_dU + cos(theta) * st_v_shear_y.v_dU;
    v_dtheta_ddU = -sin(theta) * st_v_shear_x.v_ddU + cos(theta) * st_v_shear_y.v_ddU;
    
    a_U = diag( v_U );
    a_dU = diag( v_dU );
    a_ddU = diag( v_ddU );            
    
    a_dtheta_U = diag( v_dtheta_U );
    a_dtheta_dU = diag( v_dtheta_dU );
    a_dtheta_ddU = diag( v_dtheta_ddU );
    
    % LHS
    a_LHS = a_U * a_D2m - c * a_D2m - a_ddU - k^2 * a_U + c * k^2 * a_I;
    a_LHS_FS = a_U.^2 * a_Dm - 2 * c * a_U * a_Dm + c^2 * a_Dm - a_U * a_dU + c * a_dU - ( 1 / st_p.Fr2 ) * a_I;
    a_LHS(1,:) =  a_LHS_FS(1,:);
    
    % LHS extra
    a_LHSE = -a_D2m + k^2 * a_I;
    a_LHSE_FS = -2 * a_U * a_Dm + 2 * c * a_Dm + a_dU;
    a_LHSE(1,:) = a_LHSE_FS(1,:);
    
    % RHS
    a_RHS = -a_dtheta_U * a_D2m + a_dtheta_ddU + k^2 * a_dtheta_U;
    a_RHS_FS = -2 * a_U * a_dtheta_U * a_Dm + 2 * c * a_dtheta_U * a_Dm + a_dtheta_U * a_dU + a_U * a_dtheta_dU - c * a_dtheta_dU;
    a_RHS(1,:) = a_RHS_FS(1,:);

    % Shorten
    a_LHS = a_LHS(1:end-1,1:end-1);
    a_LHSE = a_LHSE(1:end-1,1:end-1);    
    a_RHS = a_RHS(1:end-1,1:end-1);
    %v_w = v_w(1:end-1);
    
    % Full M matrix
    a_M = [ a_LHS a_LHSE*v_w; -(a_R_short*v_w)' 0 ];
    a_L = [ a_RHS * v_w; 0 ];
        
    % Solve
    warning('');  %clear warning status
    v_dy = a_M \ a_L;        
    v_dy = [ v_dy; v_dy(end) ];  %fudge
    v_dy(end-1) = 0;
    [warnMsg, warnId ] = lastwarn;
    if ( ~isempty(warnMsg) )
        error( 'Problem ill-conditioned!');
    end
    
end




function [ st_sol ] = fn__dormand_prince( t0, t1, v_y0, tol )
        
    % New setup
    tol_abs = tol;
    tol_rel = tol;
    h_facmax = 4.0;
    h_facmin = 0.05;
    h_fac = 0.925;
    h_max = 7.0;
        
    % sp=step
    maxit = 1000;
    maxit_adjust = 50;
    %h_max = 0.1;
    %h = min((t1-t0)/10, h_max );

    sgn = sign(t1-t0);
    h = sgn * 0.3;  % Empirically tested to be a reasonably start value.
    
    t = t0;
    m = length(v_y0);
    v_y = v_y0;
    
    st_sol = struct;
    st_sol.a_y = v_y;
    st_sol.a_ymid = [];
    st_sol.v_t = t;
    
    v_fy = fn__ode_f(t,v_y);
    st_sol.a_dy = v_fy;
    
    cur_itr = 0;
    feval = 0;
    adjeval = 0;
    t_calcs = 0;
%    v_h = [];
    
    while( sgn*t < sgn*t1 )
        
        if( h > h_max ), error('Step size over max'); end
        
        cur_itr = cur_itr+1;
        if(cur_itr>maxit), error('Exceeded max iterations');  end
        
        if( sgn*(t+h) >= sgn*t1-tol )
            h = t1 - t;
        end
        
        adjust_itr = 0;
        while(1)
            adjust_itr = adjust_itr + 1;
            if ( adjust_itr > maxit_adjust ), error( 'Exceeded reasonable h adjustments' ); end
            
            a_Y = [v_fy zeros(m,6)];
            for j=2:7
                v_yj = v_y + h * ( a_Y(:,1:j-1) * a_BT_Z(1:j-1,j) );
                feval = feval + 1;
                a_Y(:,j) = fn__ode_f( t + h * v_BT_c(j), v_yj );
            end
%             er = abs(h) * norm( a_Y * v_BT_b_err, 2 );
%             h_next = min( max( 0.1, min( ( ( tol / ( 2 * er ) )^0.2 ), 3 ) ) * h, h_max );
            
            % Using Hairer, Nørsett, Wanner
            v_y_diff = abs( a_Y * v_BT_b_err );
            v_y_new = v_y + h * ( a_Y * v_BT_b );
            
            v_sc = tol_abs + tol_rel * max( [ abs( v_y ), abs( v_y_new ) ], [], 2 );
            er = sqrt( (1/m) * norm( v_y_diff ./ v_sc, 2) );
            
            step_opt = (1/er)^(1/5);
            h_next = h * min( [ h_facmax, max( [ h_facmin, h_fac * step_opt ] ) ] );
            
            if ( sgn*h_next >= h_max )
                h_next = h_max - sgn*1e-7;
            end
                
            if ( er < 1 )
                break;
            end          
            
            h = h_next;
            adjeval = adjeval + 1;
%             if(er<tol)
%                 break;
%             end
%             h = min( 0.5 * h, h_next );
        end        
        v_fy = a_Y(:,end);
        v_ymid = v_y + ( h/2 ) * ( a_Y * v_BT_b_mid );
%        v_y = v_y + h * ( a_Y * v_BT_b );
        v_y = v_y_new;
        t = t + h;
        h = h_next;

%        v_h(end+1) = h;
        
        st_sol.a_y(:,end+1) = v_y;
        st_sol.a_ymid(:,end+1) = v_ymid;
        st_sol.a_dy(:,end+1) = v_fy;
        st_sol.v_t(1,end+1) = t;
        
    end    
    
    st_sol.feval = feval;
    st_sol.adjeval = adjeval;
end








% 
% 
% function [ v_dy ] = fn__ode_f( k, v_y )
%    
%     % Extract eigenvector and eigenvalue
%     v_w = v_y(1:end-1);
%     c = v_y(end);        
%         
%     % LHS
%     a_LHS = a_U * a_D2m - c * a_D2m - a_ddU - k^2 * a_U + c * k^2 * a_I;
%     a_LHS_FS = a_U^2 * a_Dm - 2 * c * a_U * a_Dm + c^2 * a_Dm - a_dU * a_U + c * a_dU - a_I;
%     a_LHS(1,:) = a_LHS_FS(1,:);    
%     
%     % LHS extra
%     a_LHSE = -a_D2m + k^2 * a_I;
%     a_LHSE_FS = -2 * a_U * a_Dm + a_dU + 2 * c * a_Dm;
%     a_LHSE(1,:) = a_LHSE_FS(1,:);
%     
%     % RHS
%     a_RHS = 2 * k * a_U - 2 * c * k * a_I;
%     a_RHS_FS = a_0I;
%     a_RHS(1,:) = a_RHS_FS(1,:);    
% 
%     % Shorten
%     a_LHS = a_LHS(1:end-1,1:end-1);
%     a_LHSE = a_LHSE(1:end-1,1:end-1);    
%     a_RHS = a_RHS(1:end-1,1:end-1);
%     v_w = v_w(1:end-1);
%     
%     % Full M matrix
%     a_M = [ a_LHS a_LHSE*v_w; -(a_R_short*v_w)' 0 ];
%     a_L = [ a_RHS * v_w; 0 ];
%         
%     % Solve
%     warning('');  %clear warning status
%     v_dy = a_M \ a_L;        
%     v_dy = [ v_dy; v_dy(end) ];  %fudge
%     v_dy(end-1) = 0;
%     [warnMsg, warnId ] = lastwarn;
%     if ( ~isempty(warnMsg) )
%         error( 'Problem ill-conditioned!');
%     end
%     
% end
% 
% 








% 
% function [ v_c_interp ] = fn__interp( st_dr, v_k )
%    
%     if( st_dr.v_k(1) > st_dr.v_k(end))
%         [~,~,P] = histcounts(-v_k,-st_dr.v_k);
%     else
%         [~,~,P] = histcounts(v_k,st_dr.v_k);
%     end
%     
%     v_c_interp = zeros(size(st_dr.v_c,1),length(P));
%       
%  
%     
%     for j=1:length(P)
%         p = P(j);
%         k1 = st_dr.v_k(p);
%         k2 = st_dr.v_k(p+1);
%         km = (k1+k2)/2;
%         dk = k2-km;
%         c1 = st_dr.v_c(:,p);
%         c2 = st_dr.v_c(:,p+1);
%         dc1 = st_dr.v_dc(:,p);
%         dc2 = st_dr.v_dc(:,p+1);
%         cm = st_dr.v_cmid(:,p);
%         s = (v_k(j)-km)/dk;
%         bf = [-1/4*s.*(2*s+3).*(s-1).^2
%               -1/4*s.*(1+s).*(s-1).^2
%               (s-1).^2.*(s+1).^2
%               -1/4*s.*(2*s-3).*(s+1).^2
%               -1/4*s.*(1-s).*(s+1).^2];
%         v_c_interp(:,j) = [c1 dk*dc1 cm c2 dk*dc2]*bf;
%     end    
%     
%     
% end
% 


end
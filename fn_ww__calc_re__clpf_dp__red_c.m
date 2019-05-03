function [ st_dp, st_err ] = fn_ww__calc_re__clpf_dp__red_c( st_Dn, k_min, k_max, st_r_shear, k0, tol, st_p )
%fn_ww__calc_re__clpf_dp__red_c: Calc PF-R-r control points (the numerical integration stage)
% 
%   [ st_dp, st_err ] = fn_ww__calc_re__clpf_dp__red_c( st_Dn, k_min, k_max, st_v_shear, k0, tol, st_p )
% 
% Performs the eigenvalue calc and numerical integration stage of the
% path-following algorithm for the reduced problem. Use
% fn_ww__calc_re__clpf_interp__sl_rad_c() to interpolate on the results.
% 
% INPUT
%   st_dp : Differentiation matrices
%   k_min, k_max : k extent required
%   st_v_shear : Struct of vector shear profile
%   k0 : Initial k point (around 2 or 3 is usually fine)
%   tol : Tolerance used for adaptive stepsize control in integrator
%   st_p : Parameter set
%   
% OUTPUT
%   st_dp : Struct containing results
%   st_err : Struct containing error handling info
%
% TAGS: CORE, SISCPFLIB
% 
% See also
%   fn_ww__calc_re__clpf_interp__sl_rad_c(),
%   fn_ww__calc_re__cl__red_c()


% Basic sanity checks
assert( isstruct( st_Dn ) );
assert( isstruct( st_r_shear ) );
assert( isstruct( st_p ) );

if ( nargout > 1 )
    st_err = struct;
    st_err.b_ok = true;
end



% After the initial asserts, we catch any errors
try

    % Process the control parameters
    if ( isfield( st_p, 'bp_pf_use_eig_mp' ) && st_p.bp_pf_use_eig_mp )
        bp_pf_use_eig_mp = st_p.bp_pf_use_eig_mp;    
        % TODO check sizes of st_Dn and st_Dn_mp match...
    else
        bp_pf_use_eig_mp = false;
    end

    % Get shear data
    [ st_v_shear, U_nd_min, U_nd_max ] = fn_ww__util__get_r_shear_data( st_Dn, st_r_shear, st_p );
    
    a_D2m = st_Dn.a_D2m;
    a_Dm = st_Dn.a_Dm;

    N = numel( st_Dn.v_z0 );


    % U matrices
    a_U = diag( st_v_shear.v_U );
    a_dU = diag( st_v_shear.v_dU );
    a_ddU = diag( st_v_shear.v_ddU );
    a_I = eye( size( a_U ) );
    a_0I = 0 * a_I;
    
    % Precompute optimisations
    a_U_D2p = a_U * a_D2m;
    a_U_Dp = a_U * a_Dm;
    a_U2_Dp = a_U^2 * a_Dm;
    a_dU_U = a_dU * a_U;
    a_U_D2p_minus_a_ddU = a_U_D2p - a_ddU;
    a_U2_Dp_minus_a_dU_U = a_U2_Dp - a_dU_U;
    a_minus2x_U_Dp_plus_dU = -2 * a_U_Dp + a_dU;


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



    % Allow initial eig calculation in mp
    if ( bp_pf_use_eig_mp )   
        [ st_p_mp ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_mp', true, 'bp_disp_update', false ) );
        [ c_CL_p, v_c_w_CL_p ] = fn_ww__calc_re__cl__red_c( st_p_mp.st_Dn_mp, [ k0 ], st_r_shear, st_p_mp );
        c_CL_p = double( c_CL_p  );
        v_c_w_CL_p = double( v_c_w_CL_p);
    else
        [ st_p_no_mp ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_disp_update', false ) );
        [ c_CL_p, v_c_w_CL_p ] = fn_ww__calc_re__cl__red_c( st_Dn, [ k0 ], st_r_shear, st_p_no_mp );
    end

    % Test for critical layer, if so throw an exception
    if ( c_CL_p <= U_nd_max )
        throw( MException( 'PFcritlayer', 'PF hit crit layer on startup' ) );
    end
    
    % Force the 2-norm to be unit
    if ( st_p.bp_pf_force_unit_evec )       
        v_c_w_CL_p = fn_ww__util__normalise_eigenvectors( v_c_w_CL_p, 2 );
    end

    % Create the eigenpair vector
    v_y0 = [ v_c_w_CL_p; c_CL_p ];


    % PF debug where we record condition number and backwards error
    if ( st_p.bp_pf_debug )
        % Cludge, this is horible.
        ret_berr = []; 
        ret_cond = [];
    end


    % Calculate the DP integral points
    [ st_sol_fw ] = fn__dormand_prince( k0, k_max, v_y0, tol );
    [ st_sol_bw ] = fn__dormand_prince( k0, k_min, v_y0, tol );


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

    v_k_PF_fw = fliplr( st_sol_fw.v_t );
    v_k_PF_bw = st_sol_bw.v_t;



    % Because the last entry of fw == first entry of bw, we need to delete one.
    % But first we take a mean, just incase there is a discrepancy.  We do not
    % however do this for the mid values, obviously.
    v_c_PF_bw_p(1) = 0.5 * ( v_c_PF_fw_p(end) + v_c_PF_bw_p(1) );
    a_w_PF_bw_p(:,1) = 0.5 * ( a_w_PF_fw_p(:,end) + a_w_PF_bw_p(:,1) );

    v_dc_PF_bw_p(1) = 0.5 * ( v_dc_PF_fw_p(end) + v_dc_PF_bw_p(1) );
    a_dw_PF_bw_p(:,1) = 0.5 * ( a_dw_PF_fw_p(:,end) + a_dw_PF_bw_p(:,1) );

    v_k_PF_bw(1) = 0.5 * ( v_k_PF_fw(end) + v_k_PF_bw(1) );

    v_c_PF_fw_p(end) = [];
    a_w_PF_fw_p(:,end) = [];

    v_dc_PF_fw_p(end) = [];
    a_dw_PF_fw_p(:,end) = [];

    v_k_PF_fw(end) = [];



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

    v_k_PF_DP = [ v_k_PF_fw v_k_PF_bw ];


    % Save everything
    st_dp = struct;
    st_dp.v_k = v_k_PF_DP;
    st_dp.v_c = v_c_PF_DP_p;
    st_dp.a_w = a_w_PF_DP_p;
    st_dp.v_dc = v_dc_PF_DP_p;
    st_dp.a_dw = a_dw_PF_DP_p;
    st_dp.v_cmid = v_cmid_PF_DP_p;
    st_dp.a_wmid = a_wmid_PF_DP_p;
    st_dp.feval_tot = feval_tot;
    st_dp.adjeval_tot = adjeval_tot;

    % Process berr and cond if required
    if ( st_p.bp_pf_debug )
        v_berr_PF_fw = fliplr( st_sol_fw.v_berr );
        v_berr_PF_bw = st_sol_bw.v_berr;

        v_cond_PF_fw = fliplr( st_sol_fw.v_cond );
        v_cond_PF_bw = st_sol_bw.v_cond;

%         v_berr_PF_bw(1) = 0.5 * ( v_berr_PF_fw(end) + v_berr_PF_bw(1) );
%         v_cond_PF_bw(1) = 0.5 * ( v_cond_PF_fw(end) + v_cond_PF_bw(1) );
% 
%         v_berr_PF_fw(end) = [];
%         v_cond_PF_fw(end) = [];

        v_berr_PF_fw(end+1) = 0.5 * ( v_berr_PF_fw(end) + v_berr_PF_bw(1) );
        v_cond_PF_fw(end+1) = 0.5 * ( v_cond_PF_fw(end) + v_cond_PF_bw(1) );

        v_berr_PF_DP = [ v_berr_PF_fw v_berr_PF_bw ];
        v_cond_PF_DP = [ v_cond_PF_fw v_cond_PF_bw ];

        st_dp.v_berr = v_berr_PF_DP;
        st_dp.v_cond = v_cond_PF_DP;
        
    end


catch ME
    
%     ME
%     ME.stack(1).name
%     ME.stack(1).line
%     ME.message 
    
    % Create a seperate function somewhere to standardise this
    if ( nargout > 1 )
        st_err.b_ok = false;
        st_dp = [];
        return;
    else        
        rethrow(ME);
    end
    
end




if ( nargout > 1 )
    st_err = struct;
    st_err.b_ok = true;   
end


% st_dp.t_ev = t_ev;
% st_dp.t_dp = t_dp;




% function [ v_y ] = fn__CL_solve( k )
%     
%     % Precomputations
%     a_A0_precomp = (  a_U * a_D2m - a_ddU  );
%     a_A2 = 0 * a_I;
% 
%     % FS
%     a_A0_FS = ( a_U^2 * a_Dm - a_dU * a_U - a_I );
%     a_A1_FS = ( -2 * a_U * a_Dm + a_dU );
%     a_A2_FS = a_Dm;
%     
%     
%     % Main equations
%     a_A0 = ( a_A0_precomp - k^2 * a_U );
%     a_A1 = -a_D2m + k^2 * a_I;
% 
%     a_A0(1,:) = a_A0_FS(1,:);
%     a_A1(1,:) = a_A1_FS(1,:);
%     a_A2(1,:) = a_A2_FS(1,:);
% 
%     a_A0(end,:) = 0 * a_A0(end,:);
%     a_A1(end,:) = 0 * a_A1(end,:);
%     a_A2(end,:) = 0 * a_A2(end,:);
%         
%     [ a_w, v_eigs ] = polyeig( a_A0, a_A1, a_A2 );
% %    v_eigs = diag( a_eig );
%     v_eigs( ~isfinite(v_eigs) ) = 0;
%     v_eigs( v_eigs < 0 ) = 0;
%     
%     % Get max eigenvalue and pull associated eigenvector
%     [ c, c_idx ] = max( v_eigs );
%     v_w = a_w( :, c_idx );
% %    v_y = [ v_w(2:end-1); mu ];   % strip off the first and last
%     v_y = [ v_w; c ];   % strip off the first and last
%     
% end




function [ v_dy ] = fn__ode_f( k, v_y )   
    
    v_w = v_y(1:end-2);    
    c = v_y(end);
            
    % LHS
    a_LHS = a_U_D2p_minus_a_ddU - c * a_D2m - k^2 * a_U + c * k^2 * a_I;
    a_LHS(1,:) = a_U2_Dp_minus_a_dU_U(1,:) - 2 * c * a_U_Dp(1,:) + c^2 * a_Dm(1,:) + c * a_dU(1,:) - ( 1 / st_p.Fr2 ) * a_I(1,:);
    
    % LHS extra
    a_LHSE = -a_D2m + k^2 * a_I;
    a_LHSE(1,:) = a_minus2x_U_Dp_plus_dU(1,:) + 2 * c * a_Dm(1,:);
    
    % RHS
    a_RHS = 2 * k * a_U - 2 * c * k * a_I;
    a_RHS(1,:) = a_0I(1,:);

    % Shorten
    a_LHS = a_LHS(1:end-1,1:end-1);
    a_LHSE = a_LHSE(1:end-1,1:end-1);    
    a_RHS = a_RHS(1:end-1,1:end-1);
    %v_w = v_w(1:end-1);
    
    % Full M matrix
    a_M = [ a_LHS a_LHSE*v_w; -(a_R_short*v_w)' 0 ];
    v_L = [ a_RHS * v_w; 0 ];
        
    % Solve
    warning('');  %clear warning status        
    v_dy = a_M \ v_L;
    
    if ( st_p.bp_pf_debug )
        ret_cond = cond( a_M );
        ret_berr = norm( v_L - a_M * v_dy, 2 ) / ( norm( a_M, 2 ) * norm( v_dy, 2 ) + norm( v_L, 2 ) );        
        %eig( a_M )
    end
    
    v_dy = [ v_dy; v_dy(end) ];  %fudge
    v_dy(end-1) = 0;
    [warnMsg, warnId ] = lastwarn;
    if ( ~isempty(warnMsg) )
        error( 'Problem ill-conditioned!');
    end
    
end






function [ st_sol ] = fn__dormand_prince( t0, t1, v_y0, tol )
        
    %fprintf( 't0=%f, t1=%f\n', t0, t1 );
    
    % New setup
    tol_abs = tol;
    tol_rel = tol;
    h_facmax = 4.0;
    h_facmin = 0.05;
    h_fac = 0.925;
    h_max = 30.0;  % TODO check this is sane
        
    % sp=step
    maxit = 1000;
    maxit_adjust = 50;
    %h_max = 0.1;
    %h = min((t1-t0)/10, h_max );

    sgn = sign(t1-t0);
    h = sgn * 0.3;  % Empirically tested to be a reasonable start value.
    
    t = t0;
    m = length(v_y0);    
    v_y = v_y0;
    
    st_sol = struct;
    st_sol.a_y = v_y;
    st_sol.a_ymid = [];
    st_sol.v_t = t;
    st_sol.v_berr = [];
    st_sol.v_cond = [];
    
    v_fy = fn__ode_f(t,v_y);
    st_sol.a_dy = v_fy;
    
    cur_itr = 0;
    feval = 0;
    adjeval = 0;
    t_calcs = 0;
%    v_h = [];
    
    while( sgn*t < sgn*t1 )
        
        if ( abs(h) > h_max )
            %warning('Step size over max');
            h = sgn*h_max;
        end
        
        cur_itr = cur_itr+1;
        %fprintf( 't=%f, cur_itr=%d, h=%f\n', t, cur_itr, h );
        if (cur_itr>maxit), error('Exceeded max iterations');  end
        
        if ( sgn*(t+h) >= sgn*t1-tol )
            h = t1 - t;
        end
        
        adjust_itr = 0;
        while(1)
            adjust_itr = adjust_itr + 1;
            if ( adjust_itr > maxit_adjust ), error( 'Exceeded reasonable h adjustments' ); end
            
            v_berr_loc = zeros( 1, 6 );
            v_cond_loc = zeros( 1, 6 );
            a_Y = [v_fy zeros(m,6)];
            for j=2:7
                v_yj = v_y + h * ( a_Y(:,1:j-1) * a_BT_Z(1:j-1,j) );
                feval = feval + 1;
                a_Y(:,j) = fn__ode_f( t + h * v_BT_c(j), v_yj );
                
                if ( st_p.bp_pf_debug )
                    v_berr_loc(j-1) = ret_berr;
                    v_cond_loc(j-1) = ret_cond;
                end
            end                                   
            
            % Using Hairer, Nørsett, Wanner
            v_y_diff = abs( a_Y * v_BT_b_err );
            v_y_new = v_y + h * ( a_Y * v_BT_b );
                                                                                          
            v_sc = tol_abs + tol_rel * max( [ abs( v_y ), abs( v_y_new ) ], [], 2 );
            er = sqrt( (1/m) * norm( v_y_diff ./ v_sc, 2 )^2 );
            
            step_opt = (1/er)^(1/5);
            h_next = h * min( [ h_facmax, max( [ h_facmin, h_fac * step_opt ] ) ] );
            
            if ( sgn*h_next >= h_max )
                h_next = h_max - sgn*1e-7;  % TODO, what does this mean?
            end
                
            if ( er < 1 )
                break;
            end          
            
            h = h_next;
            adjeval = adjeval + 1;

        end        
        
        v_fy = a_Y(:,end);    
        v_ymid = v_y + ( h/2 ) * ( a_Y * v_BT_b_mid );
        v_y = v_y_new;
        t = t + h;
        h = h_next;
        
        % Check whether in crit layer
        if ( v_y(end) <= U_nd_max )
            throw( MException( 'PFcritlayer', 'PF hit crit layer during numerical integration' ) );
        end        

        % Force the 2-norm to be unit
        if ( st_p.bp_pf_force_unit_evec )    
            [ v_y(1:end-1), scale_factor ] = fn_ww__util__normalise_eigenvectors( v_y(1:end-1), 2 );            
            v_fy(1:end-1) = v_fy(1:end-1) / scale_factor;            
            v_ymid(1:end-1) = fn_ww__util__normalise_eigenvectors( v_ymid(1:end-1), 2 );
        end
        
        st_sol.a_y(:,end+1) = v_y;
        st_sol.a_ymid(:,end+1) = v_ymid;
        st_sol.a_dy(:,end+1) = v_fy;
        st_sol.v_t(1,end+1) = t;
        if ( st_p.bp_pf_debug )
            st_sol.v_berr(end+1) = mean( v_berr_loc );
            st_sol.v_cond(end+1) = mean( v_cond_loc );
        end
    end    
    
    st_sol.feval = feval;
    st_sol.adjeval = adjeval;
end






% Before most optimisations...
% 
% function [ v_dy ] = fn__ode_f( k, v_y )
%    
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
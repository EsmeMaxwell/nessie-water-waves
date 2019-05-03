function [ v_c_p, ca_w_p, st_info ] = fn_ww__calc_shoot__red_c( v_k, st_r_shear, st_tol, st_p )
%fn_ww__calc_shoot__red_c: Calc reduced problem using shooting technique


% Shorthand
fn_U = st_r_shear.st_fn_shear.fn_U;
fn_dU = st_r_shear.st_fn_shear.fn_dU;
fn_ddU = st_r_shear.st_fn_shear.fn_ddU;

Nk = numel( v_k );


% % Tolerance set that gets 1e-10 relative error
% tol_fzero = 1e-10;
% tol_ode45_rel = 1e-8;
% tol_ode45_abs = 1e-9;
% tol_depth = 1e-10;


% % Tolerance set that gets 1e-7 relative error
% tol_fzero = 1e-6;
% tol_ode45_rel = 5e-6;
% tol_ode45_abs = 1e-6;
% tol_depth = 1e-7;

tol_fzero = st_tol.tol_fzero;
tol_ode45_rel = st_tol.tol_ode45_rel;
tol_ode45_abs = st_tol.tol_ode45_abs;
tol_depth = st_tol.tol_depth;


st_info = struct;

% Calculate initial points
[ v_z_cc, v_zw_cc ] = fn_ww__ext__cc_weights( 32, st_p.a, st_p.b );
[ v_c_EL_cc ] = fn_ww__calc_re__apx_el_cc__red_c( v_k, v_z_cc, v_zw_cc, st_r_shear.st_fn_shear, st_p );


%
ca_w_p = cell( 1, Nk );
v_zh = zeros( 1, Nk );
for lp_k=1:Nk
        
    %lp_k
    
    k = v_k( lp_k );
     
    % Calculate relevant depth
    if ( k * st_p.h < 1 )
        zh = st_p.a;
    else
        zh = log( tol_depth ) / k;
        if ( zh < st_p.a )
            zh = st_p.a;
        end
    end
    v_zh(lp_k) = zh;
    
    % Perform the shooting
    [ c, v_t, a_y ] = fn_shoot( v_c_EL_cc(lp_k) );

    % Save eigenvector data
    ca_w_p{lp_k}.v_t = v_t;
    ca_w_p{lp_k}.a_c_w = a_y( :, 1 );
    ca_w_p{lp_k}.a_c_dw = a_y( :, 2 );
    
    % Save eigenvalue data
    v_c_p( :, lp_k ) = c;   
    
end

% semilogy( v_k, v_zh )
st_info.v_zh = v_zh;



    function [ c_sht, v_t, a_y ] = fn_shoot( c0, h )

        % Do shooting
        st_options = optimset( 'TolX', tol_fzero, 'Display','iter' );
        c_sht = fzero( @fn_shoot_eval, c0 );

        % Recalculate
        dw_surface = ( fn_dU(0) * ( fn_U(0) - c_sht ) + 1/st_p.Fr2 ) / ( fn_U(0) - c_sht )^2;
        st_opts = odeset( 'RelTol', tol_ode45_rel, 'AbsTol', tol_ode45_abs, 'InitialStep', tol_ode45_abs * 0.1 );
        [ v_t, a_y ] = ode45( @( t, v_y ) fn_ODE( t, v_y, c_sht ), [ 0 zh ], [ 1 dw_surface ], st_opts );
                
    end


    function F = fn_shoot_eval( c_itr )
        
        dw_surface = ( fn_dU(0) * ( fn_U(0) - c_itr ) + 1/st_p.Fr2 ) / ( fn_U(0) - c_itr )^2;
        st_opts = odeset( 'RelTol', tol_ode45_rel, 'AbsTol', tol_ode45_abs, 'InitialStep', tol_ode45_abs * 0.1 );
        [ v_t_itr, a_y_itr ] = ode45( @( t, v_y ) fn_ODE( t, v_y, c_itr ), [ 0 zh ], [ 1 dw_surface ], st_opts );
        
        % Evaluate difference at bottom boundary
        F = a_y_itr( end, 1 );   % "- 0"
        
    end


%     function [ value_ev, isterminal_ev, direction_ev ] = fn_ODE_event( t_ev, v_y_ev )
%         
%         value_ev = 1;
%         isterminal_ev = 0;
%         direction_ev = [];
%         
%         if ( numel( v_y_ev ) > 5 && abs( v_y_ev(end) - v_y_ev(end-1) ) < 1e-7 )
%             value_ev = 0;
%             isterminal_ev = 1;
%         end                   
%             
%     end


    function [ v_dxdt ] = fn_ODE( t, v_x, c_loc )
        v_dxdt = [ 0 1; ( fn_ddU(t) / ( fn_U(t) - c_loc ) + k^2 ) 0 ] * v_x;
    end


end
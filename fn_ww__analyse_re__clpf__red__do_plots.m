function fn_ww__analyse_re__clpf__red__do_plots(  )
%fn_ww__analyse_re__clpf__red__do_plots: Plot PF 1d example output
% 
%   fn_ww__analyse_re__clpf__red__do_plots()
% 
% Plots two figures for example output from path-following algorithm: one
% normal, and one zoomed-in. Shows control and interpolation points.
% 
% TAGS: SISCPFLIB
% 
% See also
%   fn_ww__calc_re__clpf_dp__red_c(),
%   fn_ww__calc_re__clpf_interp__sl_rad_c()


% Plot stuff
set(0,'defaulttextinterpreter','latex');

% Colour defs taken from https://se.mathworks.com/help/matlab/ref/text.html
v_col_blue = [0 0.4470 0.7410];
v_col_orange = [0.8500 0.3250 0.0980];
v_col_yellow = [0.9290 0.6940 0.1250];
v_col_purple = [0.4940 0.1840 0.5560];
v_col_green = [0.4660 0.6740 0.1880];
v_col_cyan = [0.3010 0.7450 0.9330];
v_col_red =[0.6350 0.0780 0.1840];


k_min = 0.025;
k_max = 200;
k0 = 2.0;
tol = 1e-6;


[ st_p ] = fn_ww__setup__param_std__re_cl(  );
Fr2 = 0.05;
phy_h = 20;
phy_U0 = sqrt( Fr2 * st_p.phy_g * phy_h );
alpha = 2; beta = 4*pi; gamma = 1; delta = 0.5;
[ st_fn_shear_cospwr, st_p ] = fn_ww__setup__shear_fn__nondim_cospwr( st_p, phy_U0, phy_h, alpha, beta, gamma, delta );

[ st_r_shear ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear_cospwr, st_p );


% Test parameters
Nz = 48;
Nk = 80;
v_k = fn_ww__util__create_k_vec( k_min, k_max, Nk, 3, 0 );

% Setup differentiation matrix
[ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( Nz, 1 );
[ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p );


v_c_p = zeros( 1, Nk );

% Do DP integration
[ st_dp ] = fn_ww__calc_re__clpf_dp__red_c( st_Dn, k_min, k_max, st_r_shear, k0, tol, st_p );

% Do interpolation
[ v_c_PF_interp ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp, v_k, 0 );



figure(1);
plot( v_k, v_c_PF_interp, '-', 'Color', v_col_blue, 'LineWidth', 3.0 );
hold on;
scatter( st_dp.v_k, st_dp.v_c, 200, 'o', 'MarkerEdgeColor', v_col_blue, 'LineWidth', 4.0 );
% scatter( v_k, v_c_PF_interp, 60, 'r', 'x', 'LineWidth', 1.0 );

xlim( [ 0 k_max ] );
ylim( [ min( st_dp.v_c ) - 0.1 max( st_dp.v_c ) ] );
xlabel( '$k$', 'Interpreter', 'Latex', 'FontSize', 36 );
ylabel( '$c(k)$', 'Interpreter', 'Latex', 'FontSize', 36 );

set(gca,'fontsize', 36 );
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\interp_1d', '-dpdf' );




figure(2);
plot( v_k, v_c_PF_interp, '-', 'Color', v_col_blue, 'LineWidth', 4.0 );
hold on;
scatter( v_k, v_c_PF_interp, 260, '*', 'MarkerEdgeColor', v_col_red, 'LineWidth', 3.0 );
scatter( st_dp.v_k, st_dp.v_c, 200, 'o', 'MarkerEdgeColor', v_col_blue, 'LineWidth', 3.0 );


xlim( [ 10 40 ] );
ylim( [ 1.5 2.5 ] );
xlabel( '$k$', 'Interpreter', 'Latex', 'FontSize', 36 );
ylabel( '$c(k)$', 'Interpreter', 'Latex', 'FontSize', 36 );

set(gca,'fontsize', 36 );
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\interp_1d__zoom', '-dpdf' );



end
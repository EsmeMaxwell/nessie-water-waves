function [ st_gen_resample ] = fn_ww__analyse_re__clpf__2d_polar_resample__do_plots(  )
%fn_ww__analyse_re__clpf__2d_polar_resample__do_plots: DEV



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
k_max = 250;
k0 = 3.0;

tol = 1e-6;

Nz = 48;
Nthetaresample = 10;
Nkresample = 25;

Nthetaq = 37;
Nkq = 25;

v_color_lightgrey = [ 0.65 0.65 0.65 ];
v_color_darkgrey = [ 0.15 0.15 0.15 ];


k_fig1_plot_extent = 6;
k_fig21_plot_extent = 60;

v_fig1_ticks = [ -6 -4 -2 0 2 4 6 ];
v_fig21_ticks = [ -60 -40 -20 0 20 40 60 ];

v_k_fig10_xlim = [ -8 1 ];
v_k_fig10_ylim = [ -8 1 ];
v_fig10_ticks = [ -8:1:1 ];

v_z_ticks = [ -2:4:14 ];
v_z_lim = [ -2 14 ];




[ st_p ] = fn_ww__setup__param_std__re_cl(  );
[ st_fn_shear_x, st_fn_shear_y, st_p ] = fn_ww__setup__shear_fn_2d__nondim_example( st_p );

[ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( Nz, 1 );
[ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p );

[ st_v_shear_x ] = fn_ww__setup__shear_fn_to_vec( st_Dn, st_fn_shear_x, st_p );
[ st_v_shear_y ] = fn_ww__setup__shear_fn_to_vec( st_Dn, st_fn_shear_y, st_p );





[ st_gen_resample ] = fn_ww__calc_re__clpf_dp__gen_2d_polar_resample_c( st_Dn, k_min, k_max, Nkresample, Nthetaresample, st_v_shear_x, st_v_shear_y, k0, tol, st_p );


% Create full grid polar query points
v_theta_q = fn_ww__util__create_theta_vec( Nthetaq );
v_k_q = fn_ww__util__create_k_vec( k_min, k_max, Nkq, 3, 1e-4 );
[ v_theta_q_rs, v_k_q_rs, v_kx_q_rs, v_ky_q_rs ] = fn_ww__util__create_grid_ordered_pair_polar( v_k_q, v_theta_q, true );

% Create a single query point where we want to demo stage 2 calc
pt_q_kx = -3.5;
pt_q_ky = -4.0;
[ pt_q_theta, pt_q_k ] = cart2pol( pt_q_kx, pt_q_ky );

% Calc for query pts
[ v_c_PF_interp ] = fn_ww__calc_re__clpf_interp__2d_resample_c( st_gen_resample, v_theta_q_rs, v_k_q_rs );
%[ pt_q_c_interp ] = fn_ww__calc_re__clpf_interp__2d_resample_c( st_gen_resample, pt_q_theta, pt_q_k );
[ pt_q_c_interp , pt_c_radial_interp_A, pt_c_radial_interp_B, pt_theta1, pt_theta2 ] = fn_ww__dev_re__clpf_interp__2d_resample_c__debuginf( st_gen_resample, pt_q_theta, pt_q_k );

% Find pt interp location on adjoining radial slices
[ pt_q_adj1_kx, pt_q_adj1_ky ] = pol2cart( pt_theta1, pt_q_k );
[ pt_q_adj2_kx, pt_q_adj2_ky ] = pol2cart( pt_theta2, pt_q_k );


% Convert initial angular calc material into Cartesian coords for plotting
[ v_kx_theta_init_dp_angular, v_ky_theta_init_dp_angular ] = pol2cart( st_gen_resample.st_dp_angular.v_theta, st_gen_resample.st_dp_angular.k * ones( size( st_gen_resample.st_dp_angular.v_theta ) ) );

% Convert initial angular k0 at resample theta angles
[ v_kx_theta_init_interp_angular, v_ky_theta_init_interp_angular ] = pol2cart( st_gen_resample.v_theta, k0 * ones( size( st_gen_resample.v_theta ) ) );

% Convert the radial dp points into Cartesian
ca_v_kx_dp_radial = cell( 1, st_gen_resample.Ntheta);
ca_v_ky_dp_radial = cell( 1, st_gen_resample.Ntheta);
for lp_theta=1:st_gen_resample.Ntheta
    
    theta = st_gen_resample.v_theta( lp_theta );
    [ ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta} ] = pol2cart( theta * ones( size( st_gen_resample.ca_st_dp_radial{lp_theta}.v_k ) ), st_gen_resample.ca_st_dp_radial{lp_theta}.v_k );
    
end

% Convert k,theta to polar ordered pair list
[ v_theta_resample_rs, v_k_resample_rs, v_kx_resample_rs, v_ky_resample_rs ] = fn_ww__util__create_grid_ordered_pair_polar( st_gen_resample.v_k, st_gen_resample.v_theta, true );

v_c_resampled_radial = reshape( st_gen_resample.a_c_slice_radial_resampled.', 1, [] );



% Circle for fig1, fig2
Nfig1circ = 60;
v_fig1_cir_param = linspace( 0, 2*pi, Nfig1circ );
v_fig1_cir_kx = k0 * cos( v_fig1_cir_param );
v_fig1_cir_ky = k0 * sin( v_fig1_cir_param );
[ v_c_fig1circ_PF_interp ] = fn_ww__calc_re__clpf_interp__sl_ang_c( st_gen_resample.st_dp_angular, v_fig1_cir_param, 0 );

% Circle segment for fig10, fig11
Nfig10circ = 30;
v_fig10_cir_param = linspace( pt_theta1, pt_theta2, Nfig10circ );
v_fig10_cir_kx = pt_q_k * cos( v_fig10_cir_param );
v_fig10_cir_ky = pt_q_k * sin( v_fig10_cir_param );
v_fig11_proj = linspace( pt_c_radial_interp_A, pt_c_radial_interp_B, Nfig10circ );














figure(1);
%scatter( v_kx_resample_rs, v_ky_resample_rs, 260, 'r', 'x', 'LineWidth', 2.0 );
for lp_theta=1:st_gen_resample.Ntheta
    plot( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, ':', 'LineWidth', 3.0, 'Color', v_color_lightgrey );
    hold on;
%    scatter( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, 260, 'b', 'o', 'LineWidth', 2.0 );
end

%plot( [ v_kx_theta_init_interp_angular v_kx_theta_init_interp_angular(1) ], [ v_ky_theta_init_interp_angular v_ky_theta_init_interp_angular(1) ], ':', 'LineWidth', 2.0, 'Color', v_color_lightgrey );
plot( v_fig1_cir_kx, v_fig1_cir_ky, '-.', 'LineWidth', 3.0, 'Color', v_col_green );


scatter( v_kx_theta_init_dp_angular, v_ky_theta_init_dp_angular, 260, 'o', 'MarkerEdgeColor', v_col_blue, 'LineWidth', 2.0 );
scatter( v_kx_theta_init_interp_angular, v_ky_theta_init_interp_angular, 260, '*', 'MarkerEdgeColor', v_col_red, 'LineWidth', 2.0 );


xlim( [ -k_fig1_plot_extent k_fig1_plot_extent ] );
ylim( [ -k_fig1_plot_extent k_fig1_plot_extent ] );
xticks( v_fig1_ticks );
yticks( v_fig1_ticks );

xlabel( '$k_x$', 'Interpreter', 'Latex', 'Fontsize', 42 );
ylabel( '$k_y$', 'Interpreter', 'Latex', 'Fontsize', 42 );

grid on;

set(gca,'Fontsize', 42 );
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\interp_2d__precomp_ang_2d', '-dpdf' );





figure(2);
%scatter3( v_kx_resample_rs, v_ky_resample_rs, v_c_resampled_radial, 100, 'b', 'o', 'LineWidth', 1.5 );

for lp_theta=1:st_gen_resample.Ntheta
    plot3( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, st_gen_resample.ca_st_dp_radial{lp_theta}.v_c, ':', 'LineWidth', 3.0, 'Color', v_color_lightgrey );
    hold on;
%    scatter3( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, st_gen_resample.ca_st_dp_radial{lp_theta}.v_c, 160, v_color_darkgrey, '+', 'LineWidth', 1.5 );
end

v_c_init_angular = st_gen_resample.a_y0_ang_PF_interp(end,:);
%plot3( [ v_kx_theta_init_interp_angular v_kx_theta_init_interp_angular(1) ], [ v_ky_theta_init_interp_angular v_ky_theta_init_interp_angular(1) ], [ v_c_init_angular v_c_init_angular(1) ], ':', 'LineWidth', 2.0, 'Color', v_color_lightgrey );
plot3( v_fig1_cir_kx, v_fig1_cir_ky, v_c_fig1circ_PF_interp, '-.', 'LineWidth', 3.0, 'Color', v_col_green );
 


scatter3( v_kx_theta_init_dp_angular, v_ky_theta_init_dp_angular, st_gen_resample.st_dp_angular.v_c, 260, 'o', 'MarkerEdgeColor', v_col_blue, 'LineWidth', 2.0 );
scatter3( v_kx_theta_init_interp_angular, v_ky_theta_init_interp_angular, v_c_init_angular, 260, '*', 'MarkerEdgeColor', v_col_red, 'LineWidth', 2.0 );


xlim( [ -k_fig1_plot_extent k_fig1_plot_extent ] );
ylim( [ -k_fig1_plot_extent k_fig1_plot_extent ] );
zlim( v_z_lim );
xticks( v_fig1_ticks );
yticks( v_fig1_ticks );
zticks( v_z_ticks );

%view( [ -25 24 ] );
view( [ -26 15 ] );

xlabel( '$k_x$', 'Interpreter', 'Latex', 'Fontsize', 42 );
ylabel( '$k_y$', 'Interpreter', 'Latex', 'Fontsize', 42 );
zlabel( '$c$', 'Interpreter', 'Latex', 'Fontsize', 42 );

grid on;

set(gca,'fontsize', 42);
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\interp_2d__precomp_ang_3d', '-dpdf' );



















figure(3);
scatter( v_kx_resample_rs, v_ky_resample_rs, 260, 'r', '*', 'LineWidth', 2.0 );
hold on;
for lp_theta=1:st_gen_resample.Ntheta
    plot( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, ':', 'LineWidth', 2.0, 'Color', v_color_lightgrey );
    scatter( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, 260, 'o', 'MarkerEdgeColor', v_col_blue, 'LineWidth', 2.0 );
end

%plot( [ v_kx_theta_init_interp_angular v_kx_theta_init_interp_angular(1) ], [ v_ky_theta_init_interp_angular v_ky_theta_init_interp_angular(1) ], ':', 'LineWidth', 2.0, 'Color', v_color_lightgrey );
plot( v_fig1_cir_kx, v_fig1_cir_ky, '-.', 'LineWidth', 2.0, 'Color', v_col_green );


%scatter( v_kx_theta_init_dp_angular, v_ky_theta_init_dp_angular, 260, 'b', 'o', 'LineWidth', 2.0 );
scatter( v_kx_theta_init_interp_angular, v_ky_theta_init_interp_angular, 260, '*', 'MarkerEdgeColor', v_col_red, 'LineWidth', 2.0 );
scatter( v_kx_theta_init_interp_angular, v_ky_theta_init_interp_angular, 260, 'o', 'MarkerEdgeColor', v_col_blue, 'LineWidth', 2.0 );


xlim( [ -k_fig1_plot_extent k_fig1_plot_extent ] );
ylim( [ -k_fig1_plot_extent k_fig1_plot_extent ] );
xticks( v_fig1_ticks );
yticks( v_fig1_ticks );

xlabel( '$k_x$', 'Interpreter', 'Latex', 'Fontsize', 42 );
ylabel( '$k_y$', 'Interpreter', 'Latex', 'Fontsize', 42 );

grid on;

set(gca,'Fontsize', 42 );
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\interp_2d__precomp_2d', '-dpdf' );





figure(4);
scatter3( v_kx_resample_rs, v_ky_resample_rs, v_c_resampled_radial, 260, '*', 'MarkerEdgeColor', v_col_red, 'LineWidth', 2.0 );
hold on;

for lp_theta=1:st_gen_resample.Ntheta
    plot3( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, st_gen_resample.ca_st_dp_radial{lp_theta}.v_c, ':', 'LineWidth', 3.0, 'Color', v_color_lightgrey );
    scatter3( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, st_gen_resample.ca_st_dp_radial{lp_theta}.v_c, 260, 'o', 'MarkerEdgeColor', v_col_blue, 'LineWidth', 2.0 );
end

v_c_init_angular = st_gen_resample.a_y0_ang_PF_interp(end,:);
%plot3( [ v_kx_theta_init_interp_angular v_kx_theta_init_interp_angular(1) ], [ v_ky_theta_init_interp_angular v_ky_theta_init_interp_angular(1) ], [ v_c_init_angular v_c_init_angular(1) ], ':', 'LineWidth', 2.0, 'Color', v_color_lightgrey );
plot3( v_fig1_cir_kx, v_fig1_cir_ky, v_c_fig1circ_PF_interp, '-.', 'LineWidth', 3.0, 'Color', v_col_green );
 


%scatter3( v_kx_theta_init_dp_angular, v_ky_theta_init_dp_angular, st_gen_resample.st_dp_angular.v_c, 260, 'b', 'o', 'LineWidth', 2.0 );
scatter3( v_kx_theta_init_interp_angular, v_ky_theta_init_interp_angular, v_c_init_angular, 260, '*', 'MarkerEdgeColor', v_col_red, 'LineWidth', 2.0 );
scatter3( v_kx_theta_init_interp_angular, v_ky_theta_init_interp_angular, v_c_init_angular, 260, 'o', 'MarkerEdgeColor', v_col_blue, 'LineWidth', 2.0 );


xlim( [ -k_fig1_plot_extent k_fig1_plot_extent ] );
ylim( [ -k_fig1_plot_extent k_fig1_plot_extent ] );
zlim( v_z_lim );
xticks( v_fig1_ticks );
yticks( v_fig1_ticks );
zticks( v_z_ticks );

%view( [ -25 24 ] );
view( [ -26 15 ] );

xlabel( '$k_x$', 'Interpreter', 'Latex', 'Fontsize', 42 );
ylabel( '$k_y$', 'Interpreter', 'Latex', 'Fontsize', 42 );
zlabel( '$c$', 'Interpreter', 'Latex', 'Fontsize', 42 );

grid on;

set(gca,'fontsize', 42);
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\interp_2d__precomp_3d', '-dpdf' );

















% fig10 fig11 are interpolation stage
figure(10);
scatter( v_kx_resample_rs, v_ky_resample_rs, 260, '*', 'MarkerEdgeColor', v_col_red, 'LineWidth', 2.0 );
hold on;
scatter( v_kx_resample_rs, v_ky_resample_rs, 260, 'o', 'MarkerEdgeColor', v_col_blue, 'LineWidth', 2.0 );
for lp_theta=1:st_gen_resample.Ntheta
    plot( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, ':', 'LineWidth', 3.0, 'Color', v_color_lightgrey );
end


plot( v_fig10_cir_kx, v_fig10_cir_ky, '--', 'LineWidth', 3.0, 'Color', v_col_purple );

scatter( pt_q_adj1_kx, pt_q_adj1_ky, 260, 'x', 'MarkerEdgeColor', v_col_purple, 'LineWidth', 1.5 );
scatter( pt_q_adj2_kx, pt_q_adj2_ky, 260, 'x', 'MarkerEdgeColor', v_col_purple, 'LineWidth', 1.5 );

scatter( pt_q_kx, pt_q_ky, 260, '*', 'MarkerEdgeColor', v_col_red, 'LineWidth', 1.5 );

xlim( v_k_fig10_xlim );
ylim( v_k_fig10_ylim );
xticks( v_fig10_ticks );
yticks( v_fig10_ticks );

xlabel( '$k_x$', 'Interpreter', 'Latex', 'Fontsize', 42 );
ylabel( '$k_y$', 'Interpreter', 'Latex', 'Fontsize', 42 );

set(gca,'fontsize', 42);
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\interp_2d__one_pt_2d', '-dpdf' );





%pt_c_radial_interp_A
figure(11);
scatter3( v_kx_resample_rs, v_ky_resample_rs, v_c_resampled_radial, 260, '*', 'MarkerEdgeColor', v_col_red, 'LineWidth', 1.5 );
hold on;
scatter3( v_kx_resample_rs, v_ky_resample_rs, v_c_resampled_radial, 260, 'o', 'MarkerEdgeColor', v_col_blue,'LineWidth', 1.5 );
for lp_theta=1:st_gen_resample.Ntheta
    plot3( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, st_gen_resample.ca_st_dp_radial{lp_theta}.v_c, ':', 'LineWidth', 3.0, 'Color', v_color_lightgrey );
end

plot3( v_fig10_cir_kx, v_fig10_cir_ky, v_fig11_proj, '--', 'LineWidth', 3.0, 'Color', v_col_purple );

scatter3( pt_q_adj1_kx, pt_q_adj1_ky, pt_c_radial_interp_A, 260, 'x', 'MarkerEdgeColor', v_col_purple, 'LineWidth', 1.5 );
scatter3( pt_q_adj2_kx, pt_q_adj2_ky, pt_c_radial_interp_B, 260, 'x', 'MarkerEdgeColor', v_col_purple, 'LineWidth', 1.5 );

scatter3( pt_q_kx, pt_q_ky, pt_q_c_interp, 260, '*', 'MarkerEdgeColor', v_col_red, 'LineWidth', 1.5 );

xlim( v_k_fig10_xlim );
ylim( v_k_fig10_ylim );
zlim( v_z_lim );
xticks( v_fig10_ticks );
yticks( v_fig10_ticks );
zticks( v_z_ticks );
view( [ -26 15 ] );

xlabel( '$k_x$', 'Interpreter', 'Latex', 'Fontsize', 42 );
ylabel( '$k_y$', 'Interpreter', 'Latex', 'Fontsize', 42 );
zlabel( '$c$', 'Interpreter', 'Latex', 'Fontsize', 42 );

set(gca,'fontsize', 42);
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\interp_2d__one_pt_3d', '-dpdf' );







return








figure(20);
scatter( v_kx_resample_rs, v_ky_resample_rs, 100, 'b', 'o', 'LineWidth', 1.5 );
hold on;
for lp_theta=1:st_gen_resample.Ntheta
    plot( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, ':', 'LineWidth', 2.0, 'Color', v_color_lightgrey );
end

scatter( v_kx_q_rs, v_ky_q_rs, 160, 'r', '*', 'LineWidth', 1.5 );

xlim( [ -k_fig21_plot_extent k_fig21_plot_extent ] );
ylim( [ -k_fig21_plot_extent k_fig21_plot_extent ] );
xticks( v_fig21_ticks );
yticks( v_fig21_ticks );

xlabel( '$k_x$', 'Interpreter', 'Latex', 'Fontsize', 42 );
ylabel( '$k_y$', 'Interpreter', 'Latex', 'Fontsize', 42 );

set(gca,'fontsize', 42);
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\interp_2d__grid_2d', '-dpdf' );




figure(21);
scatter3( v_kx_resample_rs, v_ky_resample_rs, v_c_resampled_radial, 100, 'b', 'o', 'LineWidth', 1.5 );
hold on;
for lp_theta=1:st_gen_resample.Ntheta
    plot3( ca_v_kx_dp_radial{lp_theta}, ca_v_ky_dp_radial{lp_theta}, st_gen_resample.ca_st_dp_radial{lp_theta}.v_c, ':', 'LineWidth', 2.0, 'Color', v_color_lightgrey );
end

scatter3( v_kx_q_rs, v_ky_q_rs, v_c_PF_interp, 160, 'r', '*', 'LineWidth', 1.5 );

xlim( [ -k_fig21_plot_extent k_fig21_plot_extent ] );
ylim( [ -k_fig21_plot_extent k_fig21_plot_extent ] );
zlim( v_z_lim );
xticks( v_fig21_ticks );
yticks( v_fig21_ticks );
zticks( v_z_ticks );
view( [ -23 25 ] );

xlabel( '$k_x$', 'Interpreter', 'Latex', 'Fontsize', 42 );
ylabel( '$k_y$', 'Interpreter', 'Latex', 'Fontsize', 42 );
zlabel( '$c$', 'Interpreter', 'Latex', 'Fontsize', 42 );

set(gca,'fontsize', 42);
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\interp_2d__grid_3d', '-dpdf' );



end
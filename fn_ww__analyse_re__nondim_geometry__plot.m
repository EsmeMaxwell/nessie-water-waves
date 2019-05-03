function fn_ww__analyse_re__nondim_geometry__plot(  )
%fn_ww__analyse_re__nondim_geometry__plot: Plot physical geometry for 1st order surface wave problem
% 
%   fn_ww__analyse_re__nondim_geometry__plot()
% 
% Plots geometry for surface-wave problem used in path-following paper.
% 
% TAGS: SISCPFLIB
% 
% See also
%   fn_ww__setup__shear_fn__nondim_cospwr()



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



[ st_p ] = fn_ww__setup__param_std__re_cl(  );
Fr2 = 0.05;
phy_h = 20;
phy_U0 = sqrt( Fr2 * st_p.phy_g * phy_h );
alpha = 2; beta = 4*pi; gamma = 1; delta = 0.5;

[ st_fn_shear, st_p ] = fn_ww__setup__shear_fn__nondim_cospwr( st_p, phy_U0, phy_h, alpha, beta, gamma, delta );
[ st_r_shear ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear, st_p );

v_z_phy = linspace( st_p.phy_a, 0, 500 );


h_fig = figure(1);
plot( st_fn_shear.fn_phy_U( v_z_phy ), v_z_phy, 'Color', v_col_blue, 'LineWidth', 4.0 );
hold on;

% Labels
xlabel( 'x', 'FontSize', 36 );
ylabel( 'z', 'FontSize', 36 );

% Sort x and y ticks
xlim( [ -0.5 4 ] );
ylim( [ st_p.phy_a 1 ] );
box off;
set( gca, 'YTick', [ st_p.phy_a 0 ] );
set( gca, 'YTickLabel', { '-h', '0' } );
set( gca, 'XTick', [ 0 st_r_shear.st_crit.U_phy_min st_r_shear.st_crit.U_phy_max ] );
set( gca, 'XTickLabel', { '0', '$\acute{U}_{\textrm{min}}$', '$\acute{U}_{\textrm{max}}$' } );

set( gca, 'fontsize', 36 );
set( gca, 'TickLabelInterpreter', 'latex');

% Text label
text( 0.8, -6, '$\acute{U}(z)$', 'FontSize', 36 );


% Line
line( [ st_r_shear.st_crit.U_phy_min st_r_shear.st_crit.U_phy_min ], [ st_p.phy_a 0 ], 'LineStyle', '--', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 3.0 );
line( [ st_r_shear.st_crit.U_phy_max st_r_shear.st_crit.U_phy_max ], [ st_p.phy_a 0 ], 'LineStyle', '--', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 3.0 );
line( [ 0 0 ], [ st_p.phy_a 0 ], 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1.0 );

% Waves at top
v_x_waves = linspace( -1, 5, 500 );
line( [ -1 4 ], [ 0 0 ], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5 );
line( v_x_waves, 0.5 * sin( 3 * v_x_waves * pi ), 'Color', 'k', 'LineWidth', 3.0 );

% Sea bottom
line( [ -1 4 ], [ st_p.phy_a st_p.phy_a ], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 3.0 );

hold off;

% set( h_fig, 'Color', 'none' );
% set( h_fig, 'InvertHardCopy', 'Off' );
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\nondim_geometry', '-dpdf' );






end
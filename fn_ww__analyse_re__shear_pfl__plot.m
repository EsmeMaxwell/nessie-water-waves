function fn_ww__analyse_re__shear_pfl__plot(  )
%fn_ww__analyse_re__shear_pfl__plot: DEV
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 



set(0,'defaulttextinterpreter','latex');

Fr2 = 0.01;
phy_h = 20;


[ st_p_cospwr ] = fn_ww__setup__param_std__re_cl(  );
[ st_fn_shear_cospwr, st_p_cospwr ] = fn_ww__setup__shear_fn__nondim_cospwr( st_p_cospwr, sqrt( Fr2 * st_p_ctl.phy_g * phy_h ), phy_h, 2, 4 * pi, 1 );

[ st_p_poly ] = fn_ww__setup__param_std__re_cl(  );
[ st_fn_shear_poly, st_p_poly ] = fn_ww__setup__shear_fn__nondim_columbia_poly( st_p_poly );


v_z = linspace( -1, 0, 500 );


h_fig = figure(1);
plt_cospwr = plot( st_fn_shear_cospwr.fn_U( v_z ), v_z, 'b--', 'LineWidth', 2.0 );
hold on;
plt_poly = plot( st_fn_shear_poly.fn_U( v_z ), v_z, 'm-.', 'LineWidth', 2.0 );


% Labels
xlabel( 'x', 'FontSize', 36 );
ylabel( 'z', 'FontSize', 36 );

xlim( [ -0.4 1.2 ] );
xticks( [ -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.2 ] );
ylim( [ -1 0 ] );
yticks( [ -1 -0.8 -0.6 -0.4 -0.2 0 ] );


line( [ 0 0 ], [ -1 0 ], 'LineStyle', '-', 'Color', [ 0.1 0.1 0.1 ], 'LineWidth', 1.0 );

lgnd1 = legend( [ plt_cospwr plt_poly ], ...
        sprintf( '$U_{\\textrm{T}}$' ), ...
        sprintf( '$U_{\\textrm{CR}}$' ) );

set(lgnd1,'Location','northwest');    
set(lgnd1,'Interpreter','latex');
set(lgnd1,'FontSize', 36 );

set(gca,'fontsize', 36 );
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\shear_profile__1d', '-dpdf' );


















% 
% set(0,'defaulttextinterpreter','latex');
% 
% v_z = linspace( -1, 0, 500 );
% 
% 
% [ st_CL_fn_shear_basic ] = fn_ww__setup__shear_fn( 10, { } );
% [ st_CL_fn_shear_difficult ] = fn_ww__setup__shear_fn( 10, { 'csp_a', 1.5*pi, 'csp_b', 0.5, 'csp_c', 0 } );
% 
% 
% 
% 
% 
% 
% figure(1);
% plt_A = plot( v_z, st_CL_fn_shear_basic.fn_U( v_z ), 'b-', 'LineWidth', 2.0 );
% hold on;
% plt_B = plot( v_z, st_CL_fn_shear_difficult.fn_U( v_z ), 'g--', 'LineWidth', 2.0 );
% 
% ylim( [ -0.6 0.6 ] );
% yticks( [ -0.6 -0.4 -0.2 0 0.2 0.4 0.6 ] );
% 
% xlim( [ -1 0 ] );
% xticks( [ -1 -0.8 -0.6 -0.4 -0.2 0 ] );
% 
% %grid on;
% 
% xlabel( '$z$', 'Interpreter', 'Latex', 'FontSize', 36 );
% ylabel( '$\Upsilon(z)$', 'Interpreter', 'Latex', 'FontSize', 36 );
% 
% lgnd1 = legend( [ plt_A plt_B ], ...
%         sprintf( 'Shear profile $\\Upsilon_\\textrm{A}$' ), ...
%         sprintf( 'Shear profile $\\Upsilon_\\textrm{B}$' ) );
% 
% set(lgnd1,'Location','southwest');    
% set(lgnd1,'Interpreter','latex');
% set(lgnd1,'FontSize', 26 );
% 
% set(gca,'fontsize', 26);
% set(gca,'TickLabelInterpreter', 'latex');
% hold off;
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'PaperPositionMode','auto');
% set(gcf,'PaperUnits','normalized');
% set(gcf,'PaperPosition', [0 0 1 1]);
% print( 'figures\shear_profile__1d', '-dpdf' );
% 
% 
% 
% 
% 
% % 
% % % Difficult shear profile
% % [ st_CL_fn_shear_difficult ] = fn_ww__setup__shear_fn( 10, { 'csp_a', 1.5*pi, 'csp_b', 0.5, 'csp_c', 0 } );
% % 
% % figure(2);
% % plot( v_z, st_CL_fn_shear_difficult.fn_U( v_z ), 'b-', 'LineWidth', 2.0 );
% % hold on;
% % 
% % ylim( [ -0.6 0.6 ] );
% % yticks( [ -0.6 -0.4 -0.2 0 0.2 0.4 0.6 ] );
% % 
% % xlim( [ -1 0 ] );
% % xticks( [ -1 -0.8 -0.6 -0.4 -0.2 0 ] );
% % 
% % %grid on;
% % 
% % xlabel( '$z$', 'Interpreter', 'Latex', 'FontSize', 36 );
% % ylabel( '$U(z)$', 'Interpreter', 'Latex', 'FontSize', 36 );
% % 
% % set(gca,'fontsize', 26);
% % set(gca,'TickLabelInterpreter', 'latex');
% % hold off;
% % 
% % set(gcf,'PaperOrientation','landscape');
% % set(gcf,'PaperPositionMode','auto');
% % set(gcf,'PaperUnits','normalized');
% % set(gcf,'PaperPosition', [0 0 1 1]);
% % print( 'figures\shear_profile__difficult', '-dpdf' );
% 
% 
% 
% 
% 
% 
% % Planar shear profile
% [ st_CL_fn_shear_x ] = fn_ww__setup__shear_fn( 10, { 'csp_b', 0.05 } );
% [ st_CL_fn_shear_y ] = fn_ww__setup__shear_fn( 10, { 'csp_b', 0.025 } );
% 
% figure(10);
% plt_X = plot( v_z, st_CL_fn_shear_x.fn_U( v_z ), 'b-', 'LineWidth', 2.0 );
% hold on;
% plt_Y = plot( v_z, st_CL_fn_shear_y.fn_U( v_z ), 'g--', 'LineWidth', 2.0 );
% 
% ylim( [ 0.0 0.4 ] );
% yticks( [ 0 0.1 0.2 0.3 0.4 ] );
% 
% xlim( [ -1 0 ] );
% xticks( [ -1 -0.8 -0.6 -0.4 -0.2 0 ] );
% 
% %grid on;
% 
% xlabel( '$z$', 'Interpreter', 'Latex', 'FontSize', 36 );
% ylabel( '$\Upsilon(z)$', 'Interpreter', 'Latex', 'FontSize', 36 );
% 
% lgnd2 = legend( [ plt_X plt_Y ], ...
%         sprintf( 'Shear profile $\\Upsilon_\\textrm{Cx}$' ), ...
%         sprintf( 'Shear profile $\\Upsilon_\\textrm{Cy}$' ) );
% 
% set(lgnd2,'Location','southwest');    
% set(lgnd2,'Interpreter','latex');
% set(lgnd2,'FontSize', 26 );
% 
% set(gca,'fontsize', 26);
% set(gca,'TickLabelInterpreter', 'latex');
% hold off;
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'PaperPositionMode','auto');
% set(gcf,'PaperUnits','normalized');
% set(gcf,'PaperPosition', [0 0 1 1]);
% print( 'figures\shear_profile__2d', '-dpdf' );
% 
% 

end
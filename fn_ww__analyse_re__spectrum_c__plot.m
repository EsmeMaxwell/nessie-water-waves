function fn_ww__analyse_re__spectrum_c__plot( )
%fn_ww__analyse_re__spectrum_c__plot: Plot c spectrum for Rayleigh eqn
% 
%   fn_ww__analyse_re__spectrum_c__plot()
% 
% Plots the c spectrum for the Rayleigh eqn with free-surface boundary
% condition using the test shear profile.
% 
% TAGS: SISCPFLIB
% 
% See also
%   fn_ww__calc_re__cl__red_c__fullspec()



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


Nz = 36;
Nk = 50;
k_min = 0.25;
k_max = 250;

Fr2 = 0.05;
phy_h = 20;
alpha = 2; beta = 4*pi; gamma = 1; delta = 0.5;

% Setup
[ st_p ] = fn_ww__setup__param_std__re_cl(  );
phy_U0 = sqrt( Fr2 * st_p.phy_g * phy_h );

[ st_fn_shear, st_p ] = fn_ww__setup__shear_fn__nondim_cospwr( st_p, phy_U0, phy_h, alpha, beta, gamma, delta );
[ st_r_shear ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear, st_p );
[ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( Nz, 1 );
[ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p );
[ st_v_shear ] = fn_ww__setup__shear_fn_to_vec( st_Dn, st_fn_shear, st_p );


% Setup k vector
[ v_k ] = fn_ww__util__create_k_vec( k_min, k_max, Nk, 3, 0 );

% Calculate
[ st_p ] = fn_ww__setup__merge_parameters( st_p, struct(  ) );
[ a_c_CL_pn, a_c_CL_safe ] = fn_ww__calc_re__cl__red_c__fullspec( st_Dn, v_k, st_r_shear, st_p );


% Plot
figure(1);

% We just show continuous spectrum as a block of grey
U_min = st_r_shear.st_crit.U_min;
U_max = st_r_shear.st_crit.U_max;
rectangle( 'Position', [ k_min U_min ( k_max - k_min ) ( U_max - U_min ) ], 'FaceColor', [ 0.7 0.7 0.7 ], 'EdgeColor', [ 0.7 0.7 0.7 ] );
hold on;
 %[ k_min U_min ( k_max - k_min ) ( U_max - U_min ) ]

plot( v_k, a_c_CL_pn, 'Color', v_col_blue, 'LineWidth', 3.0 );

xlabel( '$k$', 'Interpreter', 'Latex', 'FontSize', 36 );
ylabel( '$c(k)$ in spectrum', 'Interpreter', 'Latex', 'FontSize', 36 );

xlim( [ 0 k_max ] );
ylim( [ -1 max( a_c_CL_pn(1,:) ) ] );
%ylim( [ min( a_c_CL_pn(2,:) ) max( a_c_CL_pn(1,:) ) ] );


% 
% set(lgnd5,'Location','northeast');    
% set(lgnd5,'Interpreter','latex');
% set(lgnd5,'FontSize', 26 );

set(gca,'fontsize', 26);
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\spectrum_c', '-dpdf' );




% 
% % Zoomed plot
% figure(2);
% for lp_i=1:size( a_c_CL_safe, 1 )
%     scatter( v_k, a_c_CL_safe( lp_i, : ), 3, 'MarkerEdgeColor', [ 0.7 0.7 0.7 ] );
%     hold on;
% end
% plot( v_k, a_c_CL_pn, 'b', 'LineWidth', 2.0 );
% 
% xlabel( '$k$', 'Interpreter', 'Latex', 'FontSize', 36 );
% ylabel( 'Real part of $c(k)$ in spectrum', 'Interpreter', 'Latex', 'FontSize', 36 );
% 
% xlim( [ 5 10 ] );
% ylim( [ -0.1 0.5 ] );
% %yticks( [ 1e-15 1e-13 1e-11 1e-9 1e-7 1e-5 1e-3 1e-1 ] );
% 
% % 
% % set(lgnd5,'Location','northeast');    
% % set(lgnd5,'Interpreter','latex');
% % set(lgnd5,'FontSize', 26 );
% 
% set(gca,'fontsize', 26);
% set(gca,'TickLabelInterpreter', 'latex');
% hold off;
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'PaperPositionMode','auto');
% set(gcf,'PaperUnits','normalized');
% set(gcf,'PaperPosition', [0 0 1 1]);
% print( 'figures\spectrum_c__zoom', '-dpdf' );
% 








end
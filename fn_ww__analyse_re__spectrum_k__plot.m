function [ ] = fn_ww__analyse_re__spectrum_k__plot( )
%fn_ww__analyse_re__spectrum_k__plot: Plot c spectrum for Rayleigh eqn
% 
%   fn_ww__analyse_re__spectrum_k__plot()
% 
% Plots the k spectrum for the Rayleigh eqn with free-surface boundary
% condition using the test shear profile.
% 
% TAGS: SISCPFLIB
% 
% See also
%   fn_ww__calc_re__cl__red_k__fullspec(),
%   fn_ww__calc_re__cl__red_c



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


Nz = 60;
Nk = 110;
k_min = 0.25;
k_max = 250;

Fr2 = 0.05; phy_h = 20;

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

% Calculate forward to get c vector
[ v_c ] = fn_ww__calc_re__cl__red_c( st_Dn, v_k, st_r_shear, st_p );

% Do k calculation
[ a_k_CL ] = fn_ww__calc_re__cl__red_k__fullspec( st_Dn, v_c, st_v_shear, st_p );


% Plot 
figure(1);
for lp_i=1:( size( a_k_CL, 1 ) - 1 )
    plot( v_c, a_k_CL( lp_i+1, : ), 'Color', [ 0.7 0.7 0.7 ] );
    hold on;
end
plot( v_c, a_k_CL(1,:), 'Color', v_col_blue, 'LineWidth', 3.0 );

xlabel( '$c$', 'Interpreter', 'Latex', 'FontSize', 36 );
ylabel( '$\mu(c)$', 'Interpreter', 'Latex', 'FontSize', 36 );

%xlim( [ min( v_c ) max( v_c ) ] );
%ylim( [ -0.5e4 max( a_k_CL(1,:) ) ] );
xlim( [ 1.5 max( v_c ) ] );
ylim( [ -0.5e4 0.5e4 ] );

set(gca,'fontsize', 26);
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\spectrum_mu', '-dpdf' );









% 
% 
% figure(2);
% plot( v_c, sqrt(a_k_CL(1,:)), 'b', 'LineWidth', 2.0 );
% 
% xlabel( '$c$', 'Interpreter', 'Latex', 'FontSize', 36 );
% ylabel( '$k(c)$', 'Interpreter', 'Latex', 'FontSize', 36 );
% 
% xlim( [ min( v_c ) max( v_c ) ] );
% 
% set(gca,'fontsize', 26);
% set(gca,'TickLabelInterpreter', 'latex');
% hold off;
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'PaperPositionMode','auto');
% set(gcf,'PaperUnits','normalized');
% set(gcf,'PaperPosition', [0 0 1 1]);
% print( 'figures\spectrum_k__easy', '-dpdf' );








end
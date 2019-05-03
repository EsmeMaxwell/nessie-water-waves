function [ ] = fn_ww__analyse_re__adptv__pu_plot(  )
%fn_ww__analyse_re__adptv__pu_plot: Plot partition of unity intervals
% 
%   [ ] = fn_ww__analyse_re__adptv__pu_plot(  )
% 
% Plot the partition of unity intervals for the SISC paper. 
%
% TAGS: SISCPFLIB
% 
% REQUIRES
% 
% 
% See also
%



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
Fr2 = 0.05; phy_h= 20;
[ st_fn_shear, st_p ] = fn_ww__setup__shear_fn__nondim_cospwr( st_p, sqrt( Fr2 * st_p.phy_g * phy_h ), phy_h, 2, 4 * pi, 1, 0.5 );
[ st_r_shear ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear, st_p );

[ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( 64, true );
[ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p );

[ st_Dn_dst ] = fn_ww__setup__diffmtrx__WR_poldif( 160, true );
[ st_Dn_dst ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_dst, st_p );



[ v_k ] = fn_ww__util__create_k_vec( 0.025, 250, 2000, 3, 0 );

[ st_adptv ] = fn_ww__calc_re__adptv_interval_refine( v_k, st_p );

st_proc = struct;
st_proc.ip_method = 2;
st_proc.pf_tol = 1e-11;
[ v_c_p ] = fn_ww__calc_re__adptv__red_c( st_Dn, st_r_shear, st_adptv, st_proc, st_p );


% Just manually plot it, rather than writing a loop for this...
v_k_all_idxs = 1:numel(v_k);

v_c_int1 = v_c_p;
v_c_int2 = v_c_p;
v_c_int3 = v_c_p;

v_c_int1( setdiff( v_k_all_idxs, st_adptv.ca_v_k{1}.v_k_idxs ) ) = 0 * v_c_p( setdiff( v_k_all_idxs, st_adptv.ca_v_k{1}.v_k_idxs ) );
v_c_int2( setdiff( v_k_all_idxs, st_adptv.ca_v_k{2}.v_k_idxs ) ) = 0 * v_c_p( setdiff( v_k_all_idxs, st_adptv.ca_v_k{2}.v_k_idxs ) );
v_c_int3( setdiff( v_k_all_idxs, st_adptv.ca_v_k{3}.v_k_idxs ) ) = 0 * v_c_p( setdiff( v_k_all_idxs, st_adptv.ca_v_k{3}.v_k_idxs ) );

v_c_int1( st_adptv.ca_overlap{1}.v_k_idxs ) = v_c_int1( st_adptv.ca_overlap{1}.v_k_idxs ) .* st_adptv.ca_overlap{1}.v_wt_l;
v_c_int2( st_adptv.ca_overlap{1}.v_k_idxs ) = v_c_int2( st_adptv.ca_overlap{1}.v_k_idxs ) .* st_adptv.ca_overlap{1}.v_wt_r;

v_c_int2( st_adptv.ca_overlap{2}.v_k_idxs ) = v_c_int2( st_adptv.ca_overlap{2}.v_k_idxs ) .* st_adptv.ca_overlap{2}.v_wt_l;
v_c_int3( st_adptv.ca_overlap{2}.v_k_idxs ) = v_c_int3( st_adptv.ca_overlap{2}.v_k_idxs ) .* st_adptv.ca_overlap{2}.v_wt_r;

plt_int1 = plot( v_k, v_c_int1, '--', 'Color', v_col_green, 'LineWidth', 5.0 );
hold on;
plt_int2 = plot( v_k, v_c_int2, '-.', 'Color', v_col_purple, 'LineWidth', 5.0 );
plt_int3 = plot( v_k, v_c_int3, ':', 'Color', v_col_red, 'LineWidth', 5.0 );
plt_full = plot( v_k, v_c_p, '-', 'Color', v_col_blue, 'LineWidth', 2 );


% Create interval boundaries
v_int_col = [ 0.2 0.2 0.2 ];
line( [ st_adptv.v_intervals(1,1) st_adptv.v_intervals(1,1) ], [ 5 6 ], 'LineStyle', '-', 'LineWidth', 3.0, 'Color', v_int_col );
line( [ st_adptv.v_intervals(1,2) st_adptv.v_intervals(1,2) ], [ 5 6 ], 'LineStyle', '-', 'LineWidth', 3.0, 'Color', v_int_col );
text( st_adptv.v_intervals(1,1), 5.5, '$\leftarrow$', 'FontSize', 36, 'Interpreter', 'Latex', 'Color', v_int_col );
text( st_adptv.v_intervals(1,2), 5.5, '$\rightarrow$', 'FontSize', 36, 'Interpreter', 'Latex', 'Color', v_int_col, 'HorizontalAlignment', 'right' );
text( 0.5 * ( st_adptv.v_intervals(1,1) + st_adptv.v_intervals(1,2) ), 5.5, '$I^{(1)}$', 'FontSize', 36, 'Interpreter', 'Latex', 'Color', v_int_col, 'HorizontalAlignment', 'center' );

line( [ st_adptv.v_intervals(2,1) st_adptv.v_intervals(2,1) ], [ 3.5 4.5 ], 'LineStyle', '-', 'LineWidth', 3.0, 'Color', v_int_col );
line( [ st_adptv.v_intervals(2,2) st_adptv.v_intervals(2,2) ], [ 3.5 4.5 ], 'LineStyle', '-', 'LineWidth', 3.0, 'Color', v_int_col );
text( st_adptv.v_intervals(2,1), 4, '$\leftarrow$', 'FontSize', 36, 'Interpreter', 'Latex', 'Color', v_int_col );
text( st_adptv.v_intervals(2,2), 4, '$\rightarrow$', 'FontSize', 36, 'Interpreter', 'Latex', 'Color', v_int_col, 'HorizontalAlignment', 'right' );
text( 0.5 * ( st_adptv.v_intervals(2,1) + st_adptv.v_intervals(2,2) ), 4, '$I^{(2)}$', 'FontSize', 36, 'Interpreter', 'Latex', 'Color', v_int_col, 'HorizontalAlignment', 'center' );

line( [ st_adptv.v_intervals(3,1) st_adptv.v_intervals(3,1) ], [ 2 3 ], 'LineStyle', '-', 'LineWidth', 3.0, 'Color', v_int_col );
line( [ st_adptv.v_intervals(3,2) st_adptv.v_intervals(3,2) ], [ 2 3 ], 'LineStyle', '-', 'LineWidth', 3.0, 'Color', v_int_col );
text( st_adptv.v_intervals(3,1), 2.5, '$\leftarrow$', 'FontSize', 36, 'Interpreter', 'Latex', 'Color', v_int_col );
text( st_adptv.v_intervals(3,2), 2.5, '$\rightarrow$', 'FontSize', 36, 'Interpreter', 'Latex', 'Color', v_int_col, 'HorizontalAlignment', 'right' );
text( 0.5 * ( st_adptv.v_intervals(3,1) + st_adptv.v_intervals(3,2) ), 2.5, '$I^{(3)}$', 'FontSize', 36, 'Interpreter', 'Latex', 'Color', v_int_col, 'HorizontalAlignment', 'center' );



xlabel( '$k$', 'Interpreter', 'Latex', 'FontSize', 36 );
ylabel( '$c(k)$', 'Interpreter', 'Latex', 'FontSize', 36 );


xlim( [ 0 v_k(end) ] );
ylim( [ 0 8 ] );

lgnd = legend( [ plt_int1, plt_int2, plt_int3, plt_full ], ...
            sprintf( '$I^{(1)}$ = [ %0.1f, %0.1f ]', st_adptv.v_intervals(1,1), st_adptv.v_intervals(1,2) ), ...
            sprintf( '$I^{(2)}$ = [ %0.1f, %0.1f ]', st_adptv.v_intervals(2,1), st_adptv.v_intervals(2,2) ), ...
            sprintf( '$I^{(3)}$ = [ %0.1f, %0.1f ]', st_adptv.v_intervals(3,1), st_adptv.v_intervals(3,2) ), ...
            'Full' );

set(lgnd,'Location','northeast');    
set(lgnd,'Interpreter','latex');
set(lgnd,'FontSize', 32 );

set(gca,'fontsize', 32 );
set(gca,'TickLabelInterpreter', 'latex');
hold off;

set(gcf,'PaperSize', [ 20 32 ] ); % for EGU poster
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
%print( 'figures\adptv_pu_intervals_plot', '-dpdf' );
print( 'figures\adptv_pu_intervals_plot_egu', '-dpdf' );



end
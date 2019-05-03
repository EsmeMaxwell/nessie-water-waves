function fn_ww__analyse_shear_prtb__opt_delta_h(  )
%fn_ww__analyse_shear_prtb__opt_delta_h: Plots for err-in-shear paper, optimal delta h
%
%   fn_ww__analyse_shear_prtb__basic_panels()
%
% Plots for the Water Waves paper. Uses results from
% fn_ww__sim_shear_prtb__do_runs().
% 
% TAGS: WWERRINSHEAR
%
% See also
%   fn_ww__sim_shear_prtb__do_runs()
% 

set(0,'defaulttextinterpreter','latex');

b_do_apx_check = false;

prctile_lower_idx = 2;
prctile_higher_idx = 18;

fill_alpha = 0.3;
%v_y_lim_prc = [ -3.5 3.5 ];
v_y_lim_prc1 = [ -0.15 0.15 ];
v_y_lim_prc2 = [ -0.3 0.3 ];
v_y_lim_mom = [ -0.2 0.2 ];
v_x_ticks = [ 0.3 3 30 300 ];

sz_vert_chng = 0.035;
sz_horiz_chng = 0.025;
pos_horiz_chng = 0.025;
pos_vert_chng_abs = 0.01;  %-0.035
pos_vert_chng_mult = -0.014;

% Fontsizes
fontsz_labels = 26;
fontsz_dhlabel= 18;
fontsz_axes = 26;
%fontsz_title = 24;
fontsz_title = 31;

% Linewidths
line_width_std = 3.5;
line_width_light = 3.0;
line_width_vlight = 2.0;

% Colours
col_blue = [ 34 98 247 ] / 255; 
col_blue_light = [ 160 184 239  ] / 255;

col_green = [ 42 94 70 ] / 255;
col_green_light = [ 184 224 206 ] / 255; 

col_magenta = [ 173 12 214 ] / 255;
col_magenta_light = [ 223 197 229 ] / 255;

col_red = [ 239 52 52 ] / 255;
col_red_light = [ 239 160 160 ] / 255;

col_nearblack = [ 0.1 0.1 0.1 ];

col_apx_nofilter = [ 0 0.9 0.1 ];  % Green
col_apx_filtered = [ 0.9 0.0 0.1 ];  % Red
col_pfx = [ 0.0 0.1 0.9 ];  % Blue

s_type = 'exp';

load( sprintf( 'data_shear_prtb_opt_deltah_%s/shear_prtb_%s_opt_deltah__param_space', s_type, s_type ), 'st_param_save' );



v_Fr2 = st_param_save.v_Fr2;
v_hs_prop = st_param_save.v_hs_prop;
v_delta_h_prop = st_param_save.v_delta_h_prop;
v_delta_U = st_param_save.v_delta_U;

NFr2 = numel( v_Fr2 );
Nhs = numel( v_hs_prop );
Ndeltah = numel( v_delta_h_prop );
NdeltaU = numel( v_delta_U );

ca_index = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l' };


lbl_x_base = 325;

v_h_pts = zeros( 1, Ndeltah );
v_dh_saved = zeros( 1, Ndeltah );
v_int_stddev = zeros( 1, Ndeltah );
v_int_stddev_over_num = zeros( 1, Ndeltah );
v_stddev_at_k_equals_pi = zeros( 1, Ndeltah );
i_Fr2=1;
i_counter = 0;
i_hs=1;
i_dU=1;
i_Fr2=1;
last_y = 100;
v_h_pts_to_plot = [ 10 20 50 75 100 125 250 500 ];
for lp_dh=1:1:Ndeltah
        
    i_counter = i_counter + 1;
    
    if ( v_hs_prop(i_hs) == 0 )
        b_hs = false;
    else
        b_hs = true;
    end

    st_param = st_param_save.ca_param_space{i_Fr2,i_hs,lp_dh,i_dU};           

    % Regen filename and load
    s_filename = sprintf( 'data_shear_prtb_opt_deltah_%s/shear_prtb_%s_opt_deltah__%d_%d_%d_%d', s_type, s_type, i_Fr2, i_hs, lp_dh, i_dU );
    %s_filename = sprintf( 'data_shear_prtb_opt_deltah_exp/shear_prtb_exp_opt_deltah__%d_%d_%d_%d', lp_Fr2, lp_hs, lp_dh, lp_dU );        
    load( s_filename );
    
    % Load param space file
    load( sprintf( 'data_shear_prtb_opt_deltah_%s/shear_prtb_%s_opt_deltah__param_space', s_type, s_type ) );

    % Get h_pts
    h_pts = st_param_save.ca_param_space{:,:,lp_dh}.h_pts;
    
    % Calculate relevant simulation parameters (apparently I only
    % saved some of this)
    cur_Fr2 = st_prtb_results.st_p.Fr2;
    cur_hs_prop = v_hs_prop( i_hs );
    cur_hs = st_prtb_results.st_p.h * v_hs_prop( i_hs );
    cur_delta_h_prop = v_delta_h_prop( lp_dh );
    cur_delta_h = st_prtb_results.st_p.h * v_delta_h_prop( lp_dh );
    cur_delta_U = v_delta_U( i_dU );
    cur_err_prop = ( st_prtb_results.Nitr - st_prtb_results.itr_done_pf ) / st_prtb_results.Nitr;

%         s_legend = sprintf( '$\\tilde{h}_s = %.2f$, $\\Delta \\tilde{h} = %.2f$', cur_hs, cur_delta_h );

    % Calculate integral and cost for later
    v_h_pts(lp_dh) = h_pts;
    v_dh_saved(lp_dh) = cur_delta_h;
    v_int_stddev(lp_dh) = trapz( st_prtb_results.v_k, st_prtb_results.std_err_olscb_PFX );
    v_int_stddev_over_num(lp_dh) = v_int_stddev(lp_dh) / h_pts;
   
    % We don't actually have a Pi value in our k array but we have one
    % close enough v_k(52) = 3.19
	v_stddev_at_k_equals_pi(lp_dh) = st_prtb_results.std_err_olscb_PFX(52);
    
       
%    if ( abs( last_y - st_prtb_results.std_err_olscb_PFX(end) ) > 0.00825 )
    if ( ismember( h_pts, v_h_pts_to_plot ) && h_pts == round( 1 / cur_delta_h_prop ) )
        
        fprintf( '%d, %d\n', h_pts, round( 1 / cur_delta_h_prop ) );        
        
        % Plot stddev
        h_fig_stddev = figure(1);
        semilogx( st_prtb_results.v_k, st_prtb_results.std_err_olscb_PFX, '-', 'Color', col_blue, 'LineWidth', line_width_std );
        hold on;

        % Calc text offset
        y_offset = 0;
%         if ( 10 == h_pts )
%             y_offset = 0.00215;
%         elseif ( 20 == h_pts )
%             y_offset = -0.00215;
%         end        
        if ( 10 == h_pts )
            y_offset = 0.00175;
        elseif ( 20 == h_pts )
            y_offset = -0.00175;
        end   
        
        % Disp text, there is a cludge here... we don't actually have a Pi
        % value in our k array but we have one close enough v_k(52) = 3.19
        text( lbl_x_base, st_prtb_results.std_err_olscb_PFX(end) + y_offset, sprintf( '$\\Delta \\tilde{h} = %0.3f$ ($N_h=$%d), $\\sigma(\\Delta \\tilde{c} / \\tilde{c}){|}_{k=\\pi} = %0.3f$', cur_delta_h_prop, h_pts, v_stddev_at_k_equals_pi(lp_dh) ), 'Interpreter', 'Latex', 'FontSize', fontsz_dhlabel );

        last_y = st_prtb_results.std_err_olscb_PFX(end);

    end
            
end

   

xlim( [ 0.3 300 ] );
xticks( v_x_ticks );
ylim( [ 0 0.08 ] );
% yticks( v_y_ticks_stddev );
grid( 'on' );
set( gca, 'TickLabelInterpreter', 'latex' );
set( gca, 'FontSize', fontsz_axes );
xlabel( '$\tilde{k}$', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
ylabel( sprintf( 'Stddev $\\Delta \\tilde{c} / \\tilde{c}$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
hold off;


ax_stddev = gca;
ax_stddev.Position(3) = ax_stddev.Position(3) - 0.25;
ax_stddev.Position(1) = ax_stddev.Position(1) - 0.1;
ax_stddev.Position(4) = ax_stddev.Position(4) + 0.1;
ax_stddev.Position(2) = ax_stddev.Position(2) - 0.1;

%ob_title_mean = title( sprintf( '(%s)', ca_index{lp_dU} ) );
% cur_title_pos = ob_title_mean.Position;
% cur_title_pos(1) = cur_title_pos(1) - 9.4;
% set( ob_title_mean, 'Position', cur_title_pos );

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperSize', [ 30 20 ] );
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);

s_output_filename = sprintf( 'figures_shear_prtb\\fig__%s_opt_deltah_stddev', s_type );

set( h_fig_stddev, 'Color', 'none' );
set( h_fig_stddev, 'InvertHardCopy', 'Off' );

print( s_output_filename, '-dpdf' );    

pause(1);
close( h_fig_stddev );  
   






% figure(2);
% plot( v_dh_saved, v_int_stddev, 'b' );

max_int_stddev = max( v_int_stddev );
v_int_stddev = v_int_stddev / max_int_stddev;

max_k_pi_stddev = max( v_stddev_at_k_equals_pi );
v_stddev_at_k_equals_pi = v_stddev_at_k_equals_pi / max_k_pi_stddev;

v_fit_idxs = find( v_h_pts <= 100 );
v_h_pts_fit = v_h_pts( v_fit_idxs );
v_int_stddev_fit = v_int_stddev( v_fit_idxs );
v_p_fit = polyfit( v_h_pts_fit, v_int_stddev_fit, 1 )
v_h_disp = 10:1:250;


h_fig_h_pts = figure(3);
subplot(2,1,1);
plot( v_h_pts, v_stddev_at_k_equals_pi, '.-', 'Color', col_blue, 'LineWidth', line_width_std, 'MarkerSize', 16 );
hold on;
% plot( v_h_disp, polyval( v_p_fit, v_h_disp ), '--', 'Color', col_red, 'LineWidth', line_width_std );
% text( 150, polyval( v_p_fit, 150 ), sprintf( '$%0.3f N_h + %0.3f \\qquad$', v_p_fit(1), v_p_fit(2) ), 'FontSize', fontsz_labels + 4 , 'Interpreter', 'Latex', 'HorizontalAlignment', 'right' );
xlim( [ 0 250 ] );
xticks( 0:50:250 );
grid( 'on' );
set( gca, 'TickLabelInterpreter', 'latex' );
set( gca, 'FontSize', fontsz_axes );
%xlabel( '$N_h$ (number of $h$ points)', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%ylabel( sprintf( '$\\gamma$ (integral of stddev $\\Delta \\tilde{c} / \\tilde{c})$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
ylabel( sprintf( 'Normalised $\\sigma_{N_h}(\\Delta \\tilde{c} / \\tilde{c}){|}_{k=\\pi}$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
hold off;






subplot(2,1,2);
plot( v_h_pts, v_int_stddev, '.-', 'Color', col_blue, 'LineWidth', line_width_std, 'MarkerSize', 16 );
hold on;
plot( v_h_disp, polyval( v_p_fit, v_h_disp ), '--', 'Color', col_red, 'LineWidth', line_width_std );
text( 150, polyval( v_p_fit, 150 ), sprintf( '$%0.3f N_h + %0.3f \\qquad$', v_p_fit(1), v_p_fit(2) ), 'FontSize', fontsz_labels + 4 , 'Interpreter', 'Latex', 'HorizontalAlignment', 'right' );
xlim( [ 0 250 ] );
xticks( 0:50:250 );
% xlim( [ 0 v_h_pts(end) ] );
% xticks( 0:250:v_h_pts(end) );
grid( 'on' );
set( gca, 'TickLabelInterpreter', 'latex' );
set( gca, 'FontSize', fontsz_axes );
xlabel( '$N_h$ (number of $h$ points)', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%ylabel( sprintf( '$\\gamma$ (integral of stddev $\\Delta \\tilde{c} / \\tilde{c})$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
ylabel( sprintf( '$\\gamma(N_h) / \\gamma_{\\max}$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
hold off;






set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperSize', [ 30 22 ] );
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);

s_output_filename = sprintf( 'figures_shear_prtb\\fig__%s_opt_deltah_h_pts', s_type );

set( h_fig_h_pts, 'Color', 'none' );
set( h_fig_h_pts, 'InvertHardCopy', 'Off' );

print( s_output_filename, '-dpdf' );    

pause(1);
close( h_fig_h_pts );  



end
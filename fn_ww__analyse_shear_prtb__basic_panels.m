function fn_ww__analyse_shear_prtb__basic_panels(  )
%fn_ww__analyse_shear_prtb__basic_panels: Plots for err-in-shear paper
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

% - move this plotting config stuff somewhere common...
% sz_vert_chng = 0.05;
% sz_horiz_chng =  0.04;
% pos_horiz_chng = 0.02;
% pos_vert_chng_abs = -0.008;
% pos_vert_chng_mult = -0.006;
sz_vert_chng = 0.035;
sz_horiz_chng = 0.025;
pos_horiz_chng = 0.025;
pos_vert_chng_abs = 0.01;  %-0.035
pos_vert_chng_mult = -0.014;

% Fontsizes
fontsz_labels = 32;
fontsz_axes = 28;
%fontsz_title = 24;
fontsz_title = 31;

% Linewidths
line_width_std = 3.5;
line_width_light = 3.0;

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

s_type = 'columbia';

load( sprintf( 'data_shear_prtb_%s/shear_prtb_%s__param_space', s_type, s_type ), 'st_param_save' );



v_Fr2 = st_param_save.v_Fr2;
v_hs_prop = st_param_save.v_hs_prop;
v_delta_h_prop = st_param_save.v_delta_h_prop;
v_delta_U = st_param_save.v_delta_U;

NFr2 = numel( v_Fr2 );
Nhs = numel( v_hs_prop );
Ndeltah = numel( v_delta_h_prop );
NdeltaU = numel( v_delta_U );

ca_index = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l' };


% 
% % Check how % of valid iterations for each run
% i_Fr2=2;
% if ( true ) 
% i_counter = 0;
% for lp_hs=1:Nhs
%     for lp_dh=1:Ndeltah
%         for lp_dU=1:NdeltaU
%                                  
%             % Regen filename and load
%             s_filename = sprintf( 'data_shear_prtb_%s/shear_prtb_%s__%d_%d_%d_%d', s_type, s_type, i_Fr2, lp_hs, lp_dh, lp_dU );
%             load( s_filename );
%             
%             
%             s_item = sprintf( '%s; Fr=%d, hs=%d, dh=%d, dU=%d', s_type, i_Fr2, lp_hs, lp_dh, lp_dU );
%             fprintf( '%s : %d,%d,%f\n', s_item, st_prtb_results.itr_done, st_prtb_results.Nitr, st_prtb_results.itr_done / st_prtb_results.Nitr );
% 
%         end
%     end
% end
% end        
% 
% return





% 
% % Percentiles-----
% i_Fr2=1;
% if ( true ) 
% i_counter = 0;
% for lp_hs=1:Nhs
%     for lp_dh=1:Ndeltah
%         for lp_dU=1:NdeltaU
% 
%             
%             if ( v_hs_prop(lp_hs) == 0 )
%                 b_hs = false;
%             else
%                 b_hs = true;
%             end
%             
%             i_counter = i_counter + 1;
%             
%             st_param = st_param_save.ca_param_space{i_Fr2,lp_hs,lp_dh,lp_dU};           
%             
%             % Regen filename and load
%             s_filename = sprintf( 'data_shear_prtb_%s/shear_prtb_%s__%d_%d_%d_%d', s_type, s_type, i_Fr2, lp_hs, lp_dh, lp_dU );
%             load( s_filename );
%             
%             % Calculate relevant simulation parameters (apparently I only
%             % saved some of this)
%             cur_Fr2 = st_prtb_results.st_p.Fr2;
%             cur_hs_prop = v_hs_prop( lp_hs );
%             cur_hs = st_prtb_results.st_p.h * v_hs_prop( lp_hs );
%             cur_delta_h_prop = v_delta_h_prop( lp_dh );
%             cur_delta_h = st_prtb_results.st_p.h * v_delta_h_prop( lp_dh );
%             cur_delta_U = v_delta_U( lp_dU );
%             cur_err_prop = ( st_prtb_results.Nitr - st_prtb_results.itr_done_pf ) / st_prtb_results.Nitr;
%                                    
%             
%             h_fig_percentile = figure(i_counter);
%             v_k_fill = [ st_prtb_results.v_k fliplr(st_prtb_results.v_k ) ];
% 
%             
% 
%             %--[ Percentile subplots ]-------------------------------------
%             ax_prc_sp1 = subplot( 4, 1, 1 );            
% %             semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_PF, 'Color', col_blue, 'LineWidth', line_width_std );
% %             hold on;
% %             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PF(prctile_lower_idx,:), '--', 'Color', col_blue, 'LineWidth', line_width_light );
% %             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PF(prctile_higher_idx,:), '--', 'Color', col_blue, 'LineWidth', line_width_light );
%             semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_PFX, 'Color', col_blue, 'LineWidth', line_width_std );
%             hold on;
%             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PFX(prctile_lower_idx,:), '--', 'Color', col_blue, 'LineWidth', line_width_light );
%             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PFX(prctile_higher_idx,:), '--', 'Color', col_blue, 'LineWidth', line_width_light );
%             if ( b_do_apx_check )
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_APX, '-', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_APX(prctile_lower_idx,:), '-.', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_APX(prctile_higher_idx,:), '-.', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
% 
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_APXF, '-', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_APXF(prctile_lower_idx,:), ':', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_APXF(prctile_higher_idx,:), ':', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
% 
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_PFX, '-', 'Color', col_pfx, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PFX(prctile_lower_idx,:), ':', 'Color', col_pfx, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PFX(prctile_higher_idx,:), ':', 'Color', col_pfx, 'LineWidth', line_width_std );
%             end
% %            v_olscb_PF_fill = [ st_prtb_results.a_prctl_mean_err_olscb_PF(prctile_higher_idx,:) fliplr( st_prtb_results.a_prctl_mean_err_olscb_PF(prctile_lower_idx,:) ) ];
%             v_olscb_PF_fill = [ st_prtb_results.a_prctl_mean_err_olscb_PFX(prctile_higher_idx,:) fliplr( st_prtb_results.a_prctl_mean_err_olscb_PFX(prctile_lower_idx,:) ) ];
%             h_fill_olscb_PF = fill( v_k_fill, v_olscb_PF_fill, col_blue_light );
%             set( h_fill_olscb_PF, 'FaceAlpha', fill_alpha ); %, 'EdgeAlpha', transparency );                        
%             xlim( [ 0.3 300 ] );
%             if ( ~b_hs )
%                 ylim( v_y_lim_prc1 );
%             else
%                 ylim( v_y_lim_prc2 );
%             end
%             xticks( v_x_ticks );
%             grid on;
%             %title( s_title_olscb_PF, 'Interpreter', 'Latex' );
% %            xlabel( 'k', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             set( ax_prc_sp1, 'TickLabelInterpreter', 'latex' );
%             set( ax_prc_sp1, 'FontSize', fontsz_axes );
%             set( ax_prc_sp1, 'xticklabel', [] );
%             %ylabel( 'OLS', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             text( 0.025, 0.8, '(OLS)', 'Units', 'Normalized', 'FontSize', fontsz_labels );
% %            ob_title = title( sprintf( '(%s) $\\tilde{h}_s / \\tilde{h}=%.2f$, $\\Delta \\tilde{h} / \\tilde{h}=%.2f$, $\\Delta \\tilde{U}=%.2f$', ca_index{i_counter}, cur_hs_prop, cur_delta_h_prop, cur_delta_U ), 'Interpreter', 'Latex', 'FontSize', fontsz_title );
%             ob_title = title( sprintf( '(%s) $\\tilde{h}_s =%.2f$, $\\Delta \\tilde{h} =%.2f$, $\\Delta \\tilde{U}=%.2f$', ca_index{i_counter}, cur_hs_prop, cur_delta_h_prop, cur_delta_U ), 'Interpreter', 'Latex', 'FontSize', fontsz_title );
%             if ( ~b_hs )
%                 set( ob_title, 'Position', [ 6.0 0.15 0 ] );
%             else
%                 set( ob_title, 'Position', [ 6.0 0.3 0 ] );
%             end
%             hold off;            
%                         
%             ax_prc_sp2 = subplot( 4, 1, 2 );            
% %             semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_lsqexp_PF, 'Color', col_magenta, 'LineWidth', line_width_std );
% %             hold on;            
% %             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_lsqexp_PF(prctile_lower_idx,:), '--', 'Color', col_magenta, 'LineWidth', line_width_light );
% %             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_lsqexp_PF(prctile_higher_idx,:), '--', 'Color', col_magenta, 'LineWidth', line_width_light );
%             semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_lsqexp_PFX, 'Color', col_magenta, 'LineWidth', line_width_std );
%             hold on;            
%             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_lsqexp_PFX(prctile_lower_idx,:), '--', 'Color', col_magenta, 'LineWidth', line_width_light );
%             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_lsqexp_PFX(prctile_higher_idx,:), '--', 'Color', col_magenta, 'LineWidth', line_width_light );
%             if ( b_do_apx_check )
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_lsqexp_APX, ':', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_lsqexp_APX(prctile_lower_idx,:), '-.', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_lsqexp_APX(prctile_higher_idx,:), '-.', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
% 
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_lsqexp_APXF, '-', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_lsqexp_APXF(prctile_lower_idx,:), ':', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_lsqexp_APXF(prctile_higher_idx,:), ':', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
% 
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_lsqexp_PFX, '-', 'Color', col_pfx, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_lsqexp_PFX(prctile_lower_idx,:), ':', 'Color', col_pfx, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_lsqexp_PFX(prctile_higher_idx,:), ':', 'Color', col_pfx, 'LineWidth', line_width_std );
%             end
% %             v_lsqexp_PF_fill = [ st_prtb_results.a_prctl_mean_err_lsqexp_PF(prctile_higher_idx,:) fliplr( st_prtb_results.a_prctl_mean_err_lsqexp_PF(prctile_lower_idx,:) ) ];
%             v_lsqexp_PF_fill = [ st_prtb_results.a_prctl_mean_err_lsqexp_PFX(prctile_higher_idx,:) fliplr( st_prtb_results.a_prctl_mean_err_lsqexp_PFX(prctile_lower_idx,:) ) ];
%             h_fill_lsqexp_PF = fill( v_k_fill, v_lsqexp_PF_fill, col_magenta_light );            
%             set( h_fill_lsqexp_PF, 'FaceAlpha', fill_alpha ); %, 'EdgeAlpha', transparency );            
%             xlim( [ 0.3 300 ] );
%             if ( ~b_hs )
%                 ylim( v_y_lim_prc1 );
%             else
%                 ylim( v_y_lim_prc2 );
%             end
%             xticks( v_x_ticks );
%             grid on;
%             %title( s_title_lsqexp_PF, 'Interpreter', 'Latex' );
% %            xlabel( 'k', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             set( ax_prc_sp2, 'TickLabelInterpreter', 'latex' );
%             set( ax_prc_sp2, 'FontSize', fontsz_axes );
%             set( ax_prc_sp2, 'xticklabel', [] );
%             %ylabel( 'Exp', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             text( 0.025, 0.8, '(EXP)', 'Units', 'Normalized', 'FontSize', fontsz_labels );
%             hold off;
%             
%             
%             ax_prc_sp3 = subplot( 4, 1, 3 );
% %             semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_tlscb_PF, 'Color', col_red, 'LineWidth', line_width_std );
% %             hold on;
% %             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_tlscb_PF(prctile_lower_idx,:), '--', 'Color', col_red, 'LineWidth', line_width_light );
% %             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_tlscb_PF(prctile_higher_idx,:), '--', 'Color', col_red, 'LineWidth', line_width_light );
%             semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_tlscb_PFX, 'Color', col_red, 'LineWidth', line_width_std );
%             hold on;
%             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_tlscb_PFX(prctile_lower_idx,:), '--', 'Color', col_red, 'LineWidth', line_width_light );
%             semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_tlscb_PFX(prctile_higher_idx,:), '--', 'Color', col_red, 'LineWidth', line_width_light );
%             if ( b_do_apx_check )
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_tlscb_APX, ':', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_tlscb_APX(prctile_lower_idx,:), '-.', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_tlscb_APX(prctile_higher_idx,:), '-.', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
% 
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_tlscb_APXF, '-', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_tlscb_APXF(prctile_lower_idx,:), ':', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_tlscb_APXF(prctile_higher_idx,:), ':', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
% 
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_tlscb_PFX, '-', 'Color', col_pfx, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_tlscb_PFX(prctile_lower_idx,:), ':', 'Color', col_pfx, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_tlscb_PFX(prctile_higher_idx,:), ':', 'Color', col_pfx, 'LineWidth', line_width_std );
%             end
% %             v_tlscb_PF_fill = [ st_prtb_results.a_prctl_mean_err_tlscb_PF(prctile_higher_idx,:) fliplr( st_prtb_results.a_prctl_mean_err_tlscb_PF(prctile_lower_idx,:) ) ];
%             v_tlscb_PF_fill = [ st_prtb_results.a_prctl_mean_err_tlscb_PFX(prctile_higher_idx,:) fliplr( st_prtb_results.a_prctl_mean_err_tlscb_PFX(prctile_lower_idx,:) ) ];
%             h_fill_tlscb_PF = fill( v_k_fill, v_tlscb_PF_fill, col_red_light );
%             set( h_fill_tlscb_PF, 'FaceAlpha', fill_alpha ); %, 'EdgeAlpha', transparency );  
%             xlim( [ 0.3 300 ] );
%             if ( ~b_hs )
%                 ylim( v_y_lim_prc1 );
%             else
%                 ylim( v_y_lim_prc2 );
%             end
%             xticks( v_x_ticks );            
%             grid on;
%             %title( s_title_tlscb_PF, 'Interpreter', 'Latex' );
%             set( ax_prc_sp3, 'TickLabelInterpreter', 'latex' );
%             set( ax_prc_sp3, 'FontSize', fontsz_axes );
%             if ( b_hs )
%                 set( ax_prc_sp3, 'xticklabel', [] );
%             else
%                 xlabel( '$\tilde{k}$', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             end
%             %ylabel( 'TLS', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             text( 0.025, 0.8, '(TLS)', 'Units', 'Normalized', 'FontSize', fontsz_labels );
%             hold off;
%             
%             
%             if ( b_hs )
%                 ax_prc_sp4 = subplot( 4, 1, 4 );
%                 
% %                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_PF_extr, 'Color', col_green, 'LineWidth', line_width_std );
% %                 hold on;
% %                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PF_extr(prctile_lower_idx,:), '--', 'Color', col_green, 'LineWidth', line_width_light );
% %                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PF_extr(prctile_higher_idx,:), '--', 'Color', col_green, 'LineWidth', line_width_light );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_PFX_extr, 'Color', col_green, 'LineWidth', line_width_std );
%                 hold on;
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PFX_extr(prctile_lower_idx,:), '--', 'Color', col_green, 'LineWidth', line_width_light );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PFX_extr(prctile_higher_idx,:), '--', 'Color', col_green, 'LineWidth', line_width_light );
%                 if ( b_do_apx_check )
%                     semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_APX_extr, ':', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
%                     semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_APX_extr(prctile_lower_idx,:), '-.', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
%                     semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_APX_extr(prctile_higher_idx,:), '-.', 'Color', col_apx_nofilter, 'LineWidth', line_width_std );
% 
%                     semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_APXF_extr, '-', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
%                     semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_APXF_extr(prctile_lower_idx,:), ':', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
%                     semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_APXF_extr(prctile_higher_idx,:), ':', 'Color', col_apx_filtered, 'LineWidth', line_width_std );
% 
%                     semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_PFX_extr, '-', 'Color', col_pfx, 'LineWidth', line_width_std );
%                     semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PFX_extr(prctile_lower_idx,:), ':', 'Color', col_pfx, 'LineWidth', line_width_std );
%                     semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PFX_extr(prctile_higher_idx,:), ':', 'Color', col_pfx, 'LineWidth', line_width_std );
%                 end
% %                 v_olscb_PF_fill = [ st_prtb_results.a_prctl_mean_err_olscb_PF_extr(prctile_higher_idx,:) fliplr( st_prtb_results.a_prctl_mean_err_olscb_PF_extr(prctile_lower_idx,:) ) ];
%                 v_olscb_PF_fill = [ st_prtb_results.a_prctl_mean_err_olscb_PFX_extr(prctile_higher_idx,:) fliplr( st_prtb_results.a_prctl_mean_err_olscb_PFX_extr(prctile_lower_idx,:) ) ];
%                 h_fill_olscb_PF_extr = fill( v_k_fill, v_olscb_PF_fill, col_green_light );
%                 set( h_fill_olscb_PF_extr, 'FaceAlpha', fill_alpha ); %, 'EdgeAlpha', transparency );                                
%                 xlim( [ 0.3 300 ] );
%                 if ( ~b_hs )
%                     ylim( v_y_lim_prc1 );
%                 else
%                     ylim( v_y_lim_prc2 );
%                 end
%                 xticks( v_x_ticks );
%                 grid on;
%                 %title( s_title_olscb_PF, 'Interpreter', 'Latex' );
%                 set( ax_prc_sp4, 'TickLabelInterpreter', 'latex' );
%                 set( ax_prc_sp4, 'FontSize', fontsz_axes );         
%                 xlabel( '$\tilde{k}$', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%                 %ylabel( 'Extr', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%                 text( 0.025, 0.8, '(OLS+Extr)', 'Units', 'Normalized', 'FontSize', fontsz_labels );
%                 hold off;              
%             end
%             
%             % $\Delta \tilde{c} / \tilde{c}$
%             
% %             if ( b_hs )               
% %                 text( 0.3,1.85, sprintf( '$\\tilde{h}_s=%.2f$, $\\Delta \\tilde{h}=%.2f$, $\\Delta \\tilde{U}=%.2f$', cur_hs, cur_delta_h, cur_delta_U ), 'Interpreter', 'Latex', 'FontSize', fontsz_title );
% % %                text( 0.1,1.125, sprintf( '$h_s=%.2f$ (of %.2f), $\\Delta h=%.2f$, $\\Delta U=%.2f$; err=$%.2f$', cur_hs, st_prtb_results.st_p.h, cur_delta_h, cur_delta_U, cur_err_prop ), 'Interpreter', 'Latex', 'FontSize', fontsz_title );
% %             else
% %                 text( 0.3,0.8, sprintf( '$\\tilde{h}_s=%.2f$, $\\Delta \\tilde{h}=%.2f$, $\\Delta \\tilde{U}=%.2f$', cur_hs, cur_delta_h, cur_delta_U ), 'Interpreter', 'Latex', 'FontSize', fontsz_title );
% % %                text( 0.1,0.8125, sprintf( '$h_s=%.2f$ (of %.2f), $\\Delta h=%.2f$, $\\Delta U=%.2f$; err=$%.2f$', cur_hs, st_prtb_results.st_p.h, cur_delta_h, cur_delta_U, cur_err_prop ), 'Interpreter', 'Latex', 'FontSize', fontsz_title );
% %             end
%             
%             % Hack.
%             ax_prc_sp1.Position = ax_prc_sp1.Position + [ pos_horiz_chng pos_vert_chng_abs+(3*pos_vert_chng_mult) sz_horiz_chng sz_vert_chng ];
%             ax_prc_sp2.Position = ax_prc_sp2.Position + [ pos_horiz_chng pos_vert_chng_abs+(2*pos_vert_chng_mult) sz_horiz_chng sz_vert_chng ];
%             ax_prc_sp3.Position = ax_prc_sp3.Position + [ pos_horiz_chng pos_vert_chng_abs+(1*pos_vert_chng_mult) sz_horiz_chng sz_vert_chng ];
%             if ( b_hs )
%                 ax_prc_sp4.Position = ax_prc_sp4.Position + [ pos_horiz_chng pos_vert_chng_abs sz_horiz_chng sz_vert_chng ];
%             end            
%             
%             if ( ~b_hs )
%                 p1 = get( ax_prc_sp1, 'position' );
%                 p2 = get( ax_prc_sp2, 'position' );
%                 p3 = get( ax_prc_sp3, 'position' );
%                 lbl_height = ( p1(2) + p1(4) ) - p3(2);
%                 ax_whole_y = axes( 'position', [ p3(1)-0.045 p3(2) p3(3) lbl_height], 'visible', 'off' );
%             else
%                 p1 = get( ax_prc_sp1, 'position' );
%                 p2 = get( ax_prc_sp2, 'position' );
%                 p3 = get( ax_prc_sp3, 'position' );
%                 p4 = get( ax_prc_sp4, 'position' );
%                 lbl_height = ( p1(2) + p1(4) ) - p4(2);
%                 ax_whole_y = axes( 'position', [ p4(1)-0.045 p4(2) p4(3) lbl_height], 'visible', 'off' );
%             end
%             h_whole_y_label = ylabel( '$\Delta \tilde{c} / \tilde{c}$', 'visible', 'on', 'FontSize', fontsz_labels, 'Interpreter', 'Latex' );
%             
%             %set(gca, 'color', 'none');
%             set(gcf,'PaperOrientation','portrait');  %Changed
%             set(gcf,'PaperSize', [ 20 24 ] );
%             set(gcf,'PaperPositionMode','auto');
%             set(gcf,'PaperUnits','normalized');
%             set(gcf,'PaperPosition', [0 0 1 1]);
%             
%             set(h_fig_percentile, 'Color', 'none' );
%             set(h_fig_percentile, 'InvertHardCopy', 'Off' );
%             
%             % Save
%             s_filename = sprintf( 'figures_shear_prtb\\fig__%s_shear_prctl__subplots__Fr%d_hs%d_dh%d_dU%d', s_type, i_Fr2, lp_hs, lp_dh, lp_dU );
%             s_filename_part = sprintf( 'fig__%s_shear_prctl__subplots__Fr%d_hs%d_dh%d_dU%d', s_type, i_Fr2, lp_hs, lp_dh, lp_dU);
%             disp( s_filename_part );
%             print( s_filename, '-dpdf' );
%             
%             % Now generate TeX caption for this
%             %ca_figures{( lp_Fr2 -1 )* Nhs + lp_hs,lp_dU} = sprintf( 'dU=%d, Fr2=%d, hs=%d; lp_Fr2 * Nhs + lp_hs=%d', lp_dU, lp_Fr2, lp_hs, lp_Fr2 * Nhs + lp_hs )           
%             fprintf( '\\caption{$Fr^2=%.2f$: $h_s=%.2f (%.2f)$, $\\Delta h=%.2f (%.2f)$, $\\Delta U=%.2f$; err=$%.2f$ }\n', cur_Fr2, cur_hs, cur_hs_prop, cur_delta_h, cur_delta_h_prop, cur_delta_U, cur_err_prop );
%             fprintf( '\n\n' );
%                                  
%             pause(1);
%             close( h_fig_percentile );
% 
%         end
%     end
% end
%                 
% end        
% return




















%--[ Do moments plots ]----------------------------------------------------

%v_y_lim_mean = [ -0.015 0.015 ];  -- as in orig paper
v_y_lim_mean = [ -0.15 0.15 ];  % Changed to match scale with stddev
v_y_ticks_mean = -0.15:0.05:0.15;

v_y_lim_std = [ 0 0.15 ];
v_y_ticks_std = 0:0.05:0.15;

% v_y_lim_skewness = [ -0.75 0.75]; -- as in orig paper
% v_y_ticks_skewness = [ -0.75 -0.5 -0.25 0 0.25 0.5 0.75 ];
% 
% v_y_lim_kurtosis = [ -0.75 0.75 ];
% v_y_ticks_kurtosis = [ -0.75 -0.5 -0.25 0 0.25 0.5 0.75 ];
% 

v_y_lim_skewness = [ -1.5 1.5 ];
v_y_ticks_skewness = [ -1.5 -1.0 -0.5 0 0.5 1.0 1.5 ];

v_y_lim_kurtosis = [ -1.5 1.5 ];
v_y_ticks_kurtosis = [ -1.5 -1.0 -0.5 0 0.5 1.0 1.5 ];


fontsz_legend = 26;

% Urgh.
ca_line_style = cell( 2,2 );
ca_line_style{1,1} = '-';
ca_line_style{1,2} = '--';
ca_line_style{2,1} = '-.';
ca_line_style{2,2} = ':';



i_Fr2=1;  % REMEMBER: for the exp and powerlaw, we need to set this to 2
i_counter = 0;
for lp_dU=1:NdeltaU
        
    i_counter = i_counter + 1;
    
    % Refresh figures
	h_fig_mean = figure(i_counter);    
    h_fig_std = figure(i_counter+1000);    
    h_fig_skewness = figure(i_counter+2000);
    h_fig_kurtosis = figure(i_counter+3000);
    
    for lp_hs=1:Nhs
        for lp_dh=1:Ndeltah
    
            % In this subloop we only add into existing figure

            if ( v_hs_prop(lp_hs) == 0 )
                b_hs = false;
            else
                b_hs = true;
            end

            st_param = st_param_save.ca_param_space{i_Fr2,lp_hs,lp_dh,lp_dU};           

            % Regen filename and load
            s_filename = sprintf( 'data_shear_prtb_%s/shear_prtb_%s__%d_%d_%d_%d', s_type, s_type, i_Fr2, lp_hs, lp_dh, lp_dU );
            load( s_filename );

            % Calculate relevant simulation parameters (apparently I only
            % saved some of this)
            cur_Fr2 = st_prtb_results.st_p.Fr2;
            cur_hs_prop = v_hs_prop( lp_hs );
            cur_hs = st_prtb_results.st_p.h * v_hs_prop( lp_hs );
            cur_delta_h_prop = v_delta_h_prop( lp_dh );
            cur_delta_h = st_prtb_results.st_p.h * v_delta_h_prop( lp_dh );
            cur_delta_U = v_delta_U( lp_dU );
            cur_err_prop = ( st_prtb_results.Nitr - st_prtb_results.itr_done_pf ) / st_prtb_results.Nitr;

            s_legend = sprintf( '$\\tilde{h}_s = %.2f$, $\\Delta \\tilde{h} = %.2f$', cur_hs, cur_delta_h );
            
            
            % Plot mean
            figure(i_counter);
            semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_PFX, ca_line_style{lp_hs,lp_dh}, 'Color', col_blue, 'LineWidth', line_width_std, 'DisplayName', s_legend );
            hold on;

            % Plot stddev
            figure(i_counter+1000);
            semilogx( st_prtb_results.v_k, st_prtb_results.std_err_olscb_PFX, ca_line_style{lp_hs,lp_dh}, 'Color', col_blue, 'LineWidth', line_width_std );
            hold on;
            
            % Plot skewness
            figure(i_counter+2000);
            semilogx( st_prtb_results.v_k, st_prtb_results.skew_err_olscb_PFX, ca_line_style{lp_hs,lp_dh}, 'Color', col_blue, 'LineWidth', line_width_std );
            hold on;

            % Plot kurtosis
            figure(i_counter+3000);
            semilogx( st_prtb_results.v_k, st_prtb_results.kurt_err_olscb_PFX, ca_line_style{lp_hs,lp_dh}, 'Color', col_blue, 'LineWidth', line_width_std );
            hold on;
            
        end        
        
    end
    
    % Prettify plots
    
    % Mean
    figure(i_counter);
    xlim( [ 0.3 300 ] );
    xticks( v_x_ticks );
    ylim( v_y_lim_mean );
    yticks( v_y_ticks_mean );
    grid( 'on' );
    set( gca, 'TickLabelInterpreter', 'latex' );
    set( gca, 'FontSize', fontsz_axes );
    xlabel( '$\tilde{k}$', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
	ylabel( sprintf( 'Mean $\\Delta \\tilde{c} / \\tilde{c}$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
    legend( 'Location', 'northwest', 'Orientation', 'vertical', 'Interpreter', 'Latex', 'FontSize', fontsz_legend );
    hold off;
    
    ob_title_mean = title( sprintf( '(%s)', ca_index{lp_dU} ) );
    cur_title_pos = ob_title_mean.Position;
    cur_title_pos(1) = cur_title_pos(1) - 9.4;
    set( ob_title_mean, 'Position', cur_title_pos );

    set(gcf,'PaperOrientation','portrait');  %Changed
    set(gcf,'PaperSize', [ 20 22 ] );
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperUnits','normalized');
    set(gcf,'PaperPosition', [0 0 1 1]);

    s_filename = sprintf( 'figures_shear_prtb\\fig__%s_shear_moments__mean__Fr%d_dU%d', s_type, i_Fr2, lp_dU );
    s_filename_part = sprintf( 'fig__%s_shear_moments__mean__Fr%d_dU%d', s_type, i_Fr2, lp_dU );
    
    set(h_fig_mean, 'Color', 'none' );
    set(h_fig_mean, 'InvertHardCopy', 'Off' );
    
    disp( s_filename_part );
    print( s_filename, '-dpdf' );    
                            
    pause(1);
    close( h_fig_mean );  
   

    % Stddev
    figure(i_counter+1000);
    xlim( [ 0.3 300 ] );
    xticks( v_x_ticks );
    ylim( v_y_lim_std );
    yticks( v_y_ticks_std );
    grid( 'on' );
    %legend( 'Location', 'northwest', 'Orientation', 'vertical', 'Interpreter', 'Latex', 'FontSize', fontsz_legend );
    set( gca, 'TickLabelInterpreter', 'latex' );
    set( gca, 'FontSize', fontsz_axes );
    xlabel( '$\tilde{k}$', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
	ylabel( sprintf( 'Stddev $\\Delta \\tilde{c} / \\tilde{c}$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
    hold off;
    
    ob_title_std = title( sprintf( '(%s)', ca_index{3+lp_dU} ) );
    cur_title_pos = ob_title_std.Position;
    cur_title_pos(1) = cur_title_pos(1) - 9.4;
    set( ob_title_std, 'Position', cur_title_pos );
    
    set(gcf,'PaperOrientation','portrait');  %Changed
    set(gcf,'PaperSize', [ 20 22 ] );
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperUnits','normalized');
    set(gcf,'PaperPosition', [0 0 1 1]);

    s_filename = sprintf( 'figures_shear_prtb\\fig__%s_shear_moments__stddev__Fr%d_dU%d', s_type, i_Fr2, lp_dU );
    s_filename_part = sprintf( 'fig__%s_shear_moments__stddev__Fr%d_dU%d', s_type, i_Fr2, lp_dU );
    
    set(h_fig_std, 'Color', 'none' );
    set(h_fig_std, 'InvertHardCopy', 'Off' );

    disp( s_filename_part );
    print( s_filename, '-dpdf' );    

    pause(1);
    close( h_fig_std );

    
    % Skewness
    figure(i_counter+2000);
    xlim( [ 0.3 300 ] );
    xticks( v_x_ticks );
    ylim( v_y_lim_skewness );
    yticks( v_y_ticks_skewness );
    grid( 'on' );
    %legend( 'Location', 'northwest', 'Orientation', 'vertical', 'Interpreter', 'Latex', 'FontSize', fontsz_legend );
    set( gca, 'TickLabelInterpreter', 'latex' );
    set( gca, 'FontSize', fontsz_axes );
    xlabel( '$\tilde{k}$', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
	ylabel( sprintf( 'Skewness $\\Delta \\tilde{c} / \\tilde{c}$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
    hold off;
    
    ob_title_skew = title( sprintf( '(%s)', ca_index{6+lp_dU} ) );
    cur_title_pos = ob_title_skew.Position;
    cur_title_pos(1) = cur_title_pos(1) - 9.4;
    set( ob_title_skew, 'Position', cur_title_pos );

    set(gcf,'PaperOrientation','portrait');  %Changed
    set(gcf,'PaperSize', [ 20 22 ] );
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperUnits','normalized');
    set(gcf,'PaperPosition', [0 0 1 1]);

    s_filename = sprintf( 'figures_shear_prtb\\fig__%s_shear_moments__skewness__Fr%d_dU%d', s_type, i_Fr2, lp_dU );
    s_filename_part = sprintf( 'fig__%s_shear_moments__skewness__Fr%d_dU%d', s_type, i_Fr2, lp_dU );
    
    set(h_fig_skewness, 'Color', 'none' );
    set(h_fig_skewness, 'InvertHardCopy', 'Off' );
    
    disp( s_filename_part );
    print( s_filename, '-dpdf' );    

    pause(1);
    close( h_fig_skewness );
   
    
    % Kurtosis
    figure(i_counter+3000);
    xlim( [ 0.3 300 ] );
    xticks( v_x_ticks );
    ylim( v_y_lim_kurtosis );
    yticks( v_y_ticks_kurtosis );
    grid( 'on' );
    %legend( 'Location', 'northwest', 'Orientation', 'vertical', 'Interpreter', 'Latex', 'FontSize', fontsz_legend );
    set( gca, 'TickLabelInterpreter', 'latex' );
    set( gca, 'FontSize', fontsz_axes );
    xlabel( '$\tilde{k}$', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
	ylabel( sprintf( 'Kurtosis $\\Delta \\tilde{c} / \\tilde{c}$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
    hold off;
    
    ob_title_kurt = title( sprintf( '(%s)', ca_index{9+lp_dU} ) );
    cur_title_pos = ob_title_kurt.Position;
    cur_title_pos(1) = cur_title_pos(1) - 9.4;
    set( ob_title_kurt, 'Position', cur_title_pos );

    set(gcf,'PaperOrientation','portrait');  %Changed
    set(gcf,'PaperSize', [ 20 22 ] );
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperUnits','normalized');
    set(gcf,'PaperPosition', [0 0 1 1]);

    s_filename = sprintf( 'figures_shear_prtb\\fig__%s_shear_moments__kurtosis__Fr%d_dU%d', s_type, i_Fr2, lp_dU );
    s_filename_part = sprintf( 'fig__%s_shear_moments__kurtosis__Fr%d_dU%d', s_type, i_Fr2, lp_dU );
    
    set(h_fig_kurtosis, 'Color', 'none' );
    set(h_fig_kurtosis, 'InvertHardCopy', 'Off' );    
    
    disp( s_filename_part );
    print( s_filename, '-dpdf' );    

    pause(1);
    close( h_fig_kurtosis );

    pause(1);
    
end


















            
%             
%             
%             
%             
%             %--[ Moments subplots ]----------------------------------------
%             
%             h_fig_moments  = figure(i_counter+1000);            
%             
%             ax_mom_sp1 = subplot( 4, 1, 1 );            
%             semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_PF, 'Color', col_blue, 'LineWidth', line_width_std );
%             hold on;
%             semilogx( st_prtb_results.v_k, st_prtb_results.std_err_olscb_PF, '--', 'Color', col_magenta, 'LineWidth', line_width_std );
%             semilogx( st_prtb_results.v_k, st_prtb_results.skew_err_olscb_PF, '-.', 'Color', col_green, 'LineWidth', line_width_std );
%             semilogx( st_prtb_results.v_k, st_prtb_results.kurt_err_olscb_PF, ':', 'Color', col_red, 'LineWidth', line_width_std );            
%             xlim( [ 0.3 30 ] );
%             ylim( v_y_lim_mom );
%             xticks( v_x_ticks );
%             grid on;
%             %xlabel( 'k', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             ylabel( sprintf( '\\textrm{OLS fit}\n $\\Delta c / c$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             hold off;            
%             set( ax_mom_sp1, 'TickLabelInterpreter', 'latex' );
%             set( ax_mom_sp1, 'FontSize', fontsz_axes );
%             set( ax_mom_sp1, 'xticklabel', [] );
%                  
%             
%             ax_mom_sp2 = subplot( 4, 1, 2 );            
%             semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_lsqexp_PF, 'Color', col_blue, 'LineWidth', line_width_std );
%             hold on;
%             semilogx( st_prtb_results.v_k, st_prtb_results.std_err_lsqexp_PF, '--', 'Color', col_magenta, 'LineWidth', line_width_std );
%             semilogx( st_prtb_results.v_k, st_prtb_results.skew_err_lsqexp_PF, '-.', 'Color', col_green, 'LineWidth', line_width_std );
% %            semilogx( st_prtb_results.v_k, st_prtb_results.kurt_err_lsqexp_PF, ':', 'Color', col_red, 'LineWidth', line_width_std );            
%             xlim( [ 0.3 30 ] );
%             ylim( v_y_lim_mom );
%             xticks( v_x_ticks );
%             grid on;
%             %xlabel( 'k', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             ylabel( sprintf( '\\textrm{Exp fit}\n $\\Delta c / c$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             hold off;            
%             set( ax_mom_sp2, 'TickLabelInterpreter', 'latex' );
%             set( ax_mom_sp2, 'FontSize', fontsz_axes );
%             set( ax_mom_sp2, 'xticklabel', [] );
%             
%             
%             ax_mom_sp3 = subplot( 4, 1, 3 );
%             semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_tlscb_PF, 'Color', col_blue, 'LineWidth', line_width_std );
%             hold on;
%             semilogx( st_prtb_results.v_k, st_prtb_results.std_err_tlscb_PF, '--', 'Color', col_magenta, 'LineWidth', line_width_std );
%             semilogx( st_prtb_results.v_k, st_prtb_results.skew_err_tlscb_PF, '-.', 'Color', col_green, 'LineWidth', line_width_std );
%             semilogx( st_prtb_results.v_k, st_prtb_results.kurt_err_tlscb_PF, ':', 'Color', col_red, 'LineWidth', line_width_std );            
%             xlim( [ 0.3 30 ] );
%             ylim( v_y_lim_mom );
%             xticks( v_x_ticks );
%             grid on;
%             if ( ~b_hs )
%                 xlabel( 'k', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             end
%             ylabel( sprintf( '\\textrm{TLS fit}\n $\\Delta c / c$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%             hold off;            
%             set( ax_mom_sp3, 'TickLabelInterpreter', 'latex' );
%             set( ax_mom_sp3, 'FontSize', fontsz_axes );
%             if ( b_hs )
%                 set( ax_mom_sp3, 'xticklabel', [] );
%             end
%             
%             
%             if ( b_hs )
%                 ax_mom_sp4 = subplot( 4, 1, 4 );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_PF_extr, 'Color', col_blue, 'LineWidth', line_width_std );
%                 hold on;
%                 semilogx( st_prtb_results.v_k, st_prtb_results.std_err_olscb_PF_extr, '--', 'Color', col_magenta, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.skew_err_olscb_PF_extr, '-.', 'Color', col_green, 'LineWidth', line_width_std );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.kurt_err_olscb_PF_extr, ':', 'Color', col_red, 'LineWidth', line_width_std );            
%                 xlim( [ 0.3 30 ] );
%                 ylim( v_y_lim_mom );
%                 xticks( v_x_ticks );
%                 grid on;
%                 xlabel( 'k', 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%                 ylabel( sprintf( '\\textrm{OLS+Extr}\n $\\Delta c / c$' ), 'Interpreter', 'Latex', 'FontSize', fontsz_labels );
%                 hold off;            
%                 set( ax_mom_sp4, 'TickLabelInterpreter', 'latex' );
%                 set( ax_mom_sp4, 'FontSize', fontsz_axes );                        
%             end
% 
%             % Hack.
%             ax_mom_sp1.Position = ax_mom_sp1.Position + [ 0 0 0 sz_vert_chng ];
%             ax_mom_sp2.Position = ax_mom_sp2.Position + [ 0 0 0 sz_vert_chng ];
%             ax_mom_sp3.Position = ax_mom_sp3.Position + [ 0 0 0 sz_vert_chng ];
%             if ( b_hs )
%                 ax_mom_sp4.Position = ax_mom_sp4.Position + [ 0 0 0 sz_vert_chng ];
%             end
%             
%             set(gcf,'PaperOrientation','landscape');
%             set(gcf,'PaperPositionMode','auto');
%             set(gcf,'PaperUnits','normalized');
%             set(gcf,'PaperPosition', [0 0 1 1]);
% 
%             % Save
%             s_filename = sprintf( 'figures_shear_prtb\\fig__exp_shear_moments__subplots__dU%d_Fr%d_hs%d', lp_dU, lp_Fr2, lp_hs );
%             s_filename_part = sprintf( 'fig__exp_shear_moments__dU%d_Fr%d_hs%d', lp_dU, lp_Fr2, lp_hs );
%             fprintf( 'Row Fr=%d,hs=%d, col %d\n', lp_Fr2, lp_hs, lp_dU );
%             disp( s_filename_part );
%             print( s_filename, '-dpdf' );
%             
%             
%             % Now generate TeX caption for this
%             %ca_figures{( lp_Fr2 -1 )* Nhs + lp_hs,lp_dU} = sprintf( 'dU=%d, Fr2=%d, hs=%d; lp_Fr2 * Nhs + lp_hs=%d', lp_dU, lp_Fr2, lp_hs, lp_Fr2 * Nhs + lp_hs )
%             fprintf( '\\caption{$\\Delta U=%.2f$, $Fr^2=%.2f$, $h_s=%.2f$}\n', v_delta_U(lp_dU), v_Fr2(lp_Fr2), v_hs(lp_hs) )
%             fprintf( '\n\n' );
%             
%             pause(2);
%             close( h_fig_moments );
          

% 
%                 semilogx( st_prtb_results.v_k, st_prtb_results.mean_err_olscb_PF_extr2, 'Color', [ 0.6 0.4 0.4 ], 'LineWidth', line_width_std );
%                 hold on;
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PF_extr2(prctile_lower_idx,:), '-.', 'Color', [ 0.6 0.4 0.4 ], 'LineWidth', line_width_light );
%                 semilogx( st_prtb_results.v_k, st_prtb_results.a_prctl_mean_err_olscb_PF_extr2(prctile_higher_idx,:), '-.', 'Color', [ 0.6 0.4 0.4 ], 'LineWidth', line_width_light );
%                 v_olscb_PF_fill = [ st_prtb_results.a_prctl_mean_err_olscb_PF_extr2(prctile_higher_idx,:) fliplr( st_prtb_results.a_prctl_mean_err_olscb_PF_extr2(prctile_lower_idx,:) ) ];
%                 h_fill_olscb_PF_extr2 = fill( v_k_fill, v_olscb_PF_fill, [ 0.9 0.7 0.7 ] );
%                 set( h_fill_olscb_PF_extr2, 'FaceAlpha', fill_alpha ); %, 'EdgeAlpha', transparency );

end
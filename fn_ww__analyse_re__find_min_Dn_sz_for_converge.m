function [ ab_converged_c, ab_converged_k, vb_converged_by_N_c, vb_converged_by_N_k ] = fn_ww__analyse_re__find_min_Dn_sz_for_converge( st_r_shear, st_p )
%fn_ww__analyse_re__find_min_Dn_sz_for_converge: Plot (incl. calculation) of Nz convergence by k
% 
%   [ ab_converged_c, ab_converged_k, vb_converged_by_N_c, vb_converged_by_N_k ] = fn_ww__analyse_re__find_min_Dn_sz_for_converge( st_fn_shear, st_p )
% 
% Determines required Nz matrix size for convergence at various k. This
% code is slightly fiddly and not overly reliable but it does the job.
% 
% TAGS: SISCPFLIB



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






Dn_sz_start = 12;
Dn_sz_end = 128;
% Dn_sz_start = 100;
% Dn_sz_end = 100;
Dn_sz_step = 1;

k_min = 0.025;
k_max = 250;

% thrsh_12 = 0.2;
% thrsh_9 = 0.25;
thrsh_12 = 0.16;
thrsh_9 = 0.16;


Nk = 500;

w_kc_max_deviation = 1e-4;


v_Dn_sz = Dn_sz_start:Dn_sz_step:Dn_sz_end;
N_Dn = numel( v_Dn_sz );

v_k = fn_ww__util__create_k_vec( k_min, k_max, Nk, 3, 0 );

% Vector of the first zero T coefficient
v_firstzero_T_c_idx = zeros( N_Dn, 1 );
v_firstzero_T_k_idx = zeros( N_Dn, 1 );

ab_converged_c = zeros( N_Dn, Nk );
ab_converged_k = zeros( N_Dn, Nk );    
vb_converged_by_N_c = zeros( N_Dn, 1 );
vb_converged_by_N_k = zeros( N_Dn, 1 );



for lp_N=1:N_Dn
    
    Ncur = v_Dn_sz( lp_N )
    
    % Create differentiation matrix, mapped z, shear vector
    [ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( v_Dn_sz( lp_N ), 1 );
    [ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p );
    %[ st_v_shear ] = fn_ww__setup__shear_fn_to_vec( st_Dn, st_fn_shear, st_p );    
    
    % Calculate k->c and then back c->k
    [ st_p ] = fn_ww__setup__merge_parameters( st_p, struct(  ) );
    [ v_results_CL_c, a_results_w_c ] = fn_ww__calc_re__cl__red_c( st_Dn, v_k, st_r_shear, st_p );
    [ v_results_CL_k, a_results_w_k ] = fn_ww__calc_re__cl__red_k( st_Dn, v_results_CL_c, st_r_shear, st_p );

    % Sanity check
    %[ v_k', v_results_CL_k', v_k' - v_results_CL_k' ]
%     figure(1);
%     plot( st_Dn.v_zm, a_results_w_c );
%     figure(2);
%     plot( st_Dn.v_zm, a_results_w_k );
    w_norm_check = norm( a_results_w_c - a_results_w_k, 1 );
    assert( w_norm_check < w_kc_max_deviation );
    
    % Calc T coefficients
    a_results_Tcoeff_c = fn_ww__ext__calc_T_ceoeffs( a_results_w_c );
    a_results_Tcoeff_k = fn_ww__ext__calc_T_ceoeffs( a_results_w_k );
    
    % Normalise all the coefficients to make sure they're realtive to abs
    % of maximum
    [ v_max_coeff_vals_c, v_max_coeff_idxs_c ] = max( abs( a_results_Tcoeff_c ), [], 1 );
    [ v_max_coeff_vals_k, v_max_coeff_idxs_k ] = max( abs( a_results_Tcoeff_k ), [], 1 );
    a_results_Tcoeff_c = a_results_Tcoeff_c ./ v_max_coeff_vals_c;
    a_results_Tcoeff_k = a_results_Tcoeff_k ./ v_max_coeff_vals_k;
    
    for lp_k=1:Nk
        
        % Calculate envelope a la chebfun's method
        v_envelope_c = zeros( 1, Ncur );
        v_envelope_k = zeros( 1, Ncur );
        for lp_en=1:Ncur            
            v_envelope_c(lp_en) = max( abs( a_results_Tcoeff_c( lp_en:end,lp_k ) ), [], 1 );
            v_envelope_k(lp_en) = max( abs( a_results_Tcoeff_k( lp_en:end,lp_k ) ), [], 1 );            
        end
        
        % If converging then it'll plateau, if so we will be able to find        
        % it on a suitable histogram
        [ v_hist_12_c, v_hist_edges_12_c ] = histcounts( log10( v_envelope_c ), 12  );
        v_hist_12_c = v_hist_12_c ./  v_Dn_sz( lp_N );
        [ v_hist_9_c, v_hist_edges_9_c ] = histcounts( log10( v_envelope_c ), 9  );
        v_hist_9_c = v_hist_9_c ./  v_Dn_sz( lp_N );
        
        [ v_hist_12_k, v_hist_edges_12_k ] = histcounts( log10( v_envelope_k ), 12  );
        v_hist_12_k = v_hist_12_k ./  v_Dn_sz( lp_N );
        [ v_hist_9_k, v_hist_edges_9_k ] = histcounts( log10( v_envelope_k ), 12  );
        v_hist_9_k = v_hist_9_k ./  v_Dn_sz( lp_N );
                
        % Test to see if there's a suitably long plateau, ignoring the
        % highest magnitude entry as it can have a shallow gradient
        if ( max( v_hist_12_c( 1:end-2 ) ) > thrsh_12 || max( v_hist_9_c( 1:end-2 ) ) > thrsh_9 )
            ab_converged_c( lp_N, lp_k ) = 1;
        end
        if ( max( v_hist_12_k( 1:end-2 ) ) > thrsh_12 || max( v_hist_9_k( 1:end-2 ) ) > thrsh_9 ) 
            ab_converged_k( lp_N, lp_k ) = 1;
        end
        
        % Generate plots...
        if ( false )     
            ab_converged_c( lp_N, lp_k )
            fn_local__do_plots();
            error()
        end
        
    end
    
    vb_converged_by_N_c = prod( ab_converged_c( lp_N, : ) );
    vb_converged_by_N_k = prod( ab_converged_k( lp_N, : ) );
    
end



% Plot min matrix size required for convergence against k
v_min_mtrx_sz_by_k_idx_c = zeros( 1, Nk );
v_min_mtrx_sz_by_k_c = zeros( 1, Nk );
for lp_k=1:Nk
    min_mtrx_idx = find( ab_converged_c(:,lp_k), 1, 'first' );
    if ( numel( min_mtrx_idx ) == 0 )
        min_mtrx_idx = -1;
        v_min_mtrx_sz_by_k_idx_c(lp_k) = min_mtrx_idx;
        v_min_mtrx_sz_by_k_c(lp_k) = 1024;        
    else
        v_min_mtrx_sz_by_k_idx_c(lp_k) = min_mtrx_idx;
        v_min_mtrx_sz_by_k_c(lp_k) = v_Dn_sz( v_min_mtrx_sz_by_k_idx_c(lp_k) );
    end
end




% Can't get rid of the weird spikes, so shall just envelope it for now.
v_min_mtrix_sz_by_k_c__cleaned = 0 * v_min_mtrx_sz_by_k_c;
for lp_k=1:numel(v_min_mtrx_sz_by_k_c)
    v_min_mtrix_sz_by_k_c__cleaned(lp_k) = max(v_min_mtrx_sz_by_k_c(1:lp_k));
end



figure(10);
% plot( v_k, v_min_mtrx_sz_by_k_c, 'b:', 'LineWidth', 1.0 );
% hold on;
plot( v_k, v_min_mtrix_sz_by_k_c__cleaned, 'Color', v_col_blue, 'LineWidth', 3.0 );

set(gca,'FontSize', 36 );
set(gca,'LineWidth', 2.0 );
set(gca,'TickLabelInterpreter', 'latex');

xlabel( '$k$', 'Interpreter', 'Latex','FontSize', 36  );
ylabel( 'Min reqd $N_z$', 'Interpreter', 'Latex','FontSize', 36  );

ylim( [ 0 Dn_sz_end ] );
yticks( 0:16:Dn_sz_end );
xlim( [ k_min k_max ] );

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\chebyshev_convergence_min_reqd_mtrx_sz_linear', '-dpdf' );






figure(11);
% semilogx( v_k, v_min_mtrx_sz_by_k_c, 'b:', 'LineWidth', 1.0 );
% hold on;
semilogx( v_k, v_min_mtrix_sz_by_k_c__cleaned, 'Color', v_col_blue, 'LineWidth', 3.0 );

set(gca,'FontSize', 36 );
set(gca,'LineWidth', 2.0 );
set(gca,'TickLabelInterpreter', 'latex');

xlabel( '$k$', 'Interpreter', 'Latex','FontSize', 36  );
ylabel( 'Min reqd $N_z$', 'Interpreter', 'Latex','FontSize', 36  );

ylim( [ 0 Dn_sz_end ] );
yticks( 0:16:Dn_sz_end );
xlim( [ k_min k_max ] );
xticks( [ 0.025 0.25 2.5 25 250  ] );


set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\chebyshev_convergence_min_reqd_mtrx_sz_semilog', '-dpdf' );





    function fn_local__do_plots( )
        
            % Test 1: k=end, Dn=74
            % Test 2: k=1, Dn=74
            % Test 3: k=1, Dn=96
        
            % Plot of the Chebyshev coefficients
            fig_semilogy1 = figure(1);         
            
            plt_semi1 = semilogy( abs( a_results_Tcoeff_c(:,lp_k) ), 'Color', v_col_blue, 'LineWidth', 2.0 );
            hold on;
            plt_semi2 = semilogy( v_envelope_c, 'Color', v_col_purple, 'LineWidth', 2.0 );
            hold off;
            lgnd_semi = legend( [ plt_semi1 plt_semi2 ], 'Chebyshev coefficients', 'Envelope' );
            
            set(gca,'FontSize', 26 );
            set(gca,'LineWidth', 2.0 );
            set(gca,'TickLabelInterpreter', 'latex');              
            
            set(lgnd_semi,'Location','northeast');    
            set(lgnd_semi,'Interpreter','latex');
            set(lgnd_semi,'FontSize', 26 );
            set(gca,'TickLabelInterpreter', 'latex');
            xlabel( '$n^{th}$ Chebyshev expansion coefficient', 'Interpreter', 'Latex','FontSize', 36 );
            ylabel( 'Relative magnitude of coefficient', 'Interpreter', 'Latex','FontSize', 36  );
            
            set(gcf,'PaperOrientation','landscape');
            set(gcf,'PaperPositionMode','auto');
            set(gcf,'PaperUnits','normalized');
            set(gcf,'PaperPosition', [0 0 1 1]);
            print( 'figures\chebyshev_coeffs_small_k_bigN', '-dpdf' );

            
            %[ a_results_Tcoeff_c(:,lp_k) log10( abs( a_results_Tcoeff_c(:,lp_k) ) ) ]
            
            
            % A first histogram plot
            figure(2);
            
            v_centers_12_c = (v_hist_edges_12_c(1:end-1) + v_hist_edges_12_c(2:end))/2;
            fig_bar2 = bar( v_centers_12_c, v_hist_12_c );         
            fig_bar2.FaceColor = 'flat';
            fig_bar2.CData( end, : ) = [ 0.8 0.8 0.8 ];
            fig_bar2.CData( end-1, : ) = [ 0.8 0.8 0.8 ];
                        
            % If we have convergence, colour that bar
            [ test12_val_c, test12_idx_c ] = max( v_hist_12_c( 1:end-2 ) );
            if ( test12_val_c > thrsh_12 )
                fig_bar2.CData( test12_idx_c, : ) = [ 0.8 0 0  ];
            end
      
            set(gca,'FontSize', 26 );
            set(gca,'LineWidth', 2.0 );
            set(gca,'TickLabelInterpreter', 'latex');              
            
            ylim( [ 0 0.35 ] );
            line( xlim, [ thrsh_12 thrsh_12 ], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2.0 );
            
            set(gca,'TickLabelInterpreter', 'latex');
            xlabel( 'Histogram bins', 'Interpreter', 'Latex','FontSize', 36 );
            ylabel( 'Proportion of entries in bin', 'Interpreter', 'Latex','FontSize', 36 );
                      
            set(gcf,'PaperOrientation','landscape');
            set(gcf,'PaperPositionMode','auto');
            set(gcf,'PaperUnits','normalized');
            set(gcf,'PaperPosition', [0 0 1 1]);
            print( 'figures\chebyshev_hist_small_k_bigN', '-dpdf' );     
            
            figure(3);
            v_centers_9_c = (v_hist_edges_9_c(1:end-1) + v_hist_edges_9_c(2:end))/2;
            fig_bar3 = bar( v_centers_9_c, v_hist_9_c );    
     
    end



end
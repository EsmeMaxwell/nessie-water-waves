function [ ] = fn_ww__analyse_re__perf_tests_c__red(  )
%fn_ww__analyse_re__perf_tests_c__red: Plot performance results for CL,DIM,PF
% 
%   fn_ww__analyse_re__perf_tests_c__red()
% 
% Plot the performance results for CL, DIM, PF precomputed by
% fn_ww__sim_re__perf_tests_c__red__set().
%
% TAGS: SISCPFLIB
% 
% REQUIRES
% 
% Files sim_perf_c_*
% 
% See also
%   fn_ww__sim_re__perf_tests_c__red__set(),
%   fn_ww__sim_re__perf_tests_c__red(),
%   fn_ww__sim_re__perf_tests_c__red__parfor()



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


s_path = 'sisc_data';

ca_s_files = cell( 21, 1 );
ca_s_files{1} = sprintf( '%s/sim_perf_c_set29__Nk100_ref__k0_025_250__cospwr', s_path );
ca_s_files{2} = sprintf( '%s/sim_perf_c_set29__Nk300_50_1000__k0_025_250__cospwr', s_path );
ca_s_files{3} = sprintf( '%s/sim_perf_c_set29__Nk1125_125_2500__k0_025_250__cospwr', s_path );
ca_s_files{4} = sprintf( '%s/sim_perf_c_set29__Nk3000_250_5000__k0_025_250__cospwr', s_path );
ca_s_files{5} = sprintf( '%s/sim_perf_c_set29__Nk6000_500_7500__k0_025_250__cospwr', s_path );
ca_s_files{6} = sprintf( '%s/sim_perf_c_set29__Nk8000_1000_10000__k0_025_250__cospwr', s_path );
ca_s_files{7} = sprintf( '%s/sim_perf_c_set29__Nk11000_1000_12000__k0_025_250__cospwr', s_path );
ca_s_files{8} = sprintf( '%s/sim_perf_c_set29__Nk13000_1000_14000__k0_025_250__cospwr', s_path );
ca_s_files{9} = sprintf( '%s/sim_perf_c_set29__Nk15000_1000_16000__k0_025_250__cospwr', s_path );
ca_s_files{10} = sprintf( '%s/sim_perf_c_set29__Nk17000_1000_18000__k0_025_250__cospwr', s_path );
ca_s_files{11} = sprintf( '%s/sim_perf_c_set29__Nk19000_1000_20000__k0_025_250__cospwr', s_path );
ca_s_files{12} = sprintf( '%s/sim_perf_c_set29__Nk21000_1000_30000__k0_025_250__cospwr', s_path );
ca_s_files{13} = sprintf( '%s/sim_perf_c_set29__Nk32500_2500_50000__k0_025_250__cospwr', s_path );
ca_s_files{14} = sprintf( '%s/sim_perf_c_set29__Nk55000_5000_100000__k0_025_250__cospwr', s_path );
ca_s_files{15} = sprintf( '%s/sim_perf_c_set29__Nk150000_10000_200000__k0_025_250__cospwr', s_path );
ca_s_files{16} = sprintf( '%s/sim_perf_c_set29__Nk250000_50000_500000__k0_025_250__cospwr', s_path );
ca_s_files{17} = sprintf( '%s/sim_perf_c_set29__Nk600000_100000_1000000__k0_025_250__cospwr', s_path );
ca_s_files{18} = sprintf( '%s/sim_perf_c_set29__Nk1100000_100000_1500000__k0_025_250__cospwr', s_path );
ca_s_files{19} = sprintf( '%s/sim_perf_c_set29__Nk1600000_100000_2000000__k0_025_250__cospwr', s_path );
ca_s_files{20} = sprintf( '%s/sim_perf_c_set29__Nk2100000_100000_2500000__k0_025_250__cospwr', s_path );
ca_s_files{21} = sprintf( '%s/sim_perf_c_set29__Nk2600000_100000_3000000__k0_025_250__cospwr', s_path );





v_t_CL_c_low = [ ];
v_t_CL_c_med = [ ];
v_t_CL_c_high = [ ];

v_t_SH_c_low = [ ];
v_t_SH_c_high = [ ];

v_t_PF_c_low = [ ];
v_t_PF_c_med = [ ];
v_t_PF_c_high = [ ];

v_t_ADPTV_c_low = [ ];
v_t_ADPTV_c_med = [ ];
v_t_ADPTV_c_high = [ ];

v_t_PFmp_c_low = [ ];
v_t_PFmp_c_med = [ ];
v_t_PFmp_c_high = [ ];

v_t_DIM_low = [ ];
v_t_DIM_med = [ ];
v_t_DIM_high = [ ];



v_all_CL_c_mp_ref = [ ];

v_all_CL_c_16 = [ ];
v_all_CL_c_low = [ ];
v_all_CL_c_med = [ ];
v_all_CL_c_high = [ ];

v_all_SH_c_low = [ ];
v_all_SH_c_high = [ ];

v_all_PF_c_low = [ ];
v_all_PF_c_med = [ ];
v_all_PF_c_high = [ ];

v_all_ADPTV_c_low = [ ];
v_all_ADPTV_c_med = [ ];
v_all_ADPTV_c_high = [ ];

v_all_PFmp_c_low = [ ];
v_all_PFmp_c_med = [ ];
v_all_PFmp_c_high = [ ];

v_all_DIM_low = [ ];
v_all_DIM_med = [ ];
v_all_DIM_high = [ ];

v_Nk = [ ];





for lp_f=1:numel(ca_s_files)

    fprintf( 'Processing file %d of %d.\n', lp_f, numel(ca_s_files) );
    
    % Load
    load( ca_s_files{lp_f} );    
    
    % Collcate index
    v_Nk = [ v_Nk st_perf_tests_c.v_Nk ];
    
    % Collate time stuff
    v_t_CL_c_low = [ v_t_CL_c_low; st_perf_tests_c.v_t_CL_c_low ];
    v_t_CL_c_med = [ v_t_CL_c_med; st_perf_tests_c.v_t_CL_c_med ];
    v_t_CL_c_high = [ v_t_CL_c_high; st_perf_tests_c.v_t_CL_c_high ];
    
    v_t_SH_c_low = [ v_t_SH_c_low; st_perf_tests_c.v_t_SH_c_low.' ];   % TODO fix in original output
    v_t_SH_c_high = [ v_t_SH_c_high; st_perf_tests_c.v_t_SH_c_high.' ];
    
    v_t_PF_c_low = [ v_t_PF_c_low; st_perf_tests_c.v_t_PF_c_low ];
    v_t_PF_c_med = [ v_t_PF_c_med; st_perf_tests_c.v_t_PF_c_med ];
    v_t_PF_c_high = [ v_t_PF_c_high; st_perf_tests_c.v_t_PF_c_high ];
    
    v_t_ADPTV_c_low = [ v_t_ADPTV_c_low; st_perf_tests_c.v_t_ADPTV_c_low ];
    v_t_ADPTV_c_med = [ v_t_ADPTV_c_med; st_perf_tests_c.v_t_ADPTV_c_med ];
    v_t_ADPTV_c_high = [ v_t_ADPTV_c_high; st_perf_tests_c.v_t_ADPTV_c_high ];

    v_t_PFmp_c_low = [ v_t_PFmp_c_low; st_perf_tests_c.v_t_PFmp_c_low ];
    v_t_PFmp_c_med = [ v_t_PFmp_c_med; st_perf_tests_c.v_t_PFmp_c_med ];
    v_t_PFmp_c_high = [ v_t_PFmp_c_high; st_perf_tests_c.v_t_PFmp_c_high ];

    v_t_DIM_low = [ v_t_DIM_low; st_perf_tests_c.v_t_DIM_low ];
    v_t_DIM_med = [ v_t_DIM_med; st_perf_tests_c.v_t_DIM_med ];
    v_t_DIM_high = [ v_t_DIM_high; st_perf_tests_c.v_t_DIM_high ];
        
    Nk_segment = numel( st_perf_tests_c.v_Nk );

    if ( st_perf_tests_c.ca_v_c_CL_mp_ref{1} > 0 )
    
        for lp_Nk=1:Nk_segment

            fprintf( 'Nk-segment %d of %d.\n', lp_Nk, Nk_segment );
            
            % This section is horrific for memory allocation but I don't
            % think it can be avoided with the split files.
            v_all_CL_c_mp_ref = [ v_all_CL_c_mp_ref, st_perf_tests_c.ca_v_c_CL_mp_ref{lp_Nk} ];
            
            v_all_CL_c_low = [ v_all_CL_c_low, st_perf_tests_c.ca_v_c_CL_low{lp_Nk} ];
            v_all_CL_c_med = [ v_all_CL_c_med, st_perf_tests_c.ca_v_c_CL_med{lp_Nk} ];
            v_all_CL_c_high = [ v_all_CL_c_high, st_perf_tests_c.ca_v_c_CL_high{lp_Nk} ];
            
            v_all_SH_c_low = [ v_all_SH_c_low, st_perf_tests_c.ca_v_c_SH_low{lp_Nk} ];
            v_all_SH_c_high = [ v_all_SH_c_high, st_perf_tests_c.ca_v_c_SH_high{lp_Nk} ];

            v_all_PF_c_low = [ v_all_PF_c_low, st_perf_tests_c.ca_v_c_PF_interp_low{lp_Nk} ];
            v_all_PF_c_med = [ v_all_PF_c_med, st_perf_tests_c.ca_v_c_PF_interp_med{lp_Nk} ];
            v_all_PF_c_high = [ v_all_PF_c_high, st_perf_tests_c.ca_v_c_PF_interp_high{lp_Nk} ];

            v_all_ADPTV_c_low = [ v_all_ADPTV_c_low, st_perf_tests_c.ca_v_c_ADPTV_low{lp_Nk} ];
            v_all_ADPTV_c_med = [ v_all_ADPTV_c_med, st_perf_tests_c.ca_v_c_ADPTV_med{lp_Nk} ];
            v_all_ADPTV_c_high = [ v_all_ADPTV_c_high, st_perf_tests_c.ca_v_c_ADPTV_high{lp_Nk} ];

            v_all_PFmp_c_low = [ v_all_PFmp_c_low, st_perf_tests_c.ca_v_c_PFmp_interp_low{lp_Nk} ];
            v_all_PFmp_c_med = [ v_all_PFmp_c_med, st_perf_tests_c.ca_v_c_PFmp_interp_med{lp_Nk} ];
            v_all_PFmp_c_high = [ v_all_PFmp_c_high, st_perf_tests_c.ca_v_c_PFmp_interp_high{lp_Nk} ];

            v_all_DIM_low = [ v_all_DIM_low, st_perf_tests_c.ca_v_tc_dimM_low{lp_Nk} ];
            v_all_DIM_med = [ v_all_DIM_med, st_perf_tests_c.ca_v_tc_dimM_med{lp_Nk} ];
            v_all_DIM_high = [ v_all_DIM_high, st_perf_tests_c.ca_v_tc_dimM_high{lp_Nk} ];
            
        end              
        
    else
        
        fprintf( 'Not using error data from this file, %d, as no CL data (this is ok).\n', lp_f );
        
    end % End error calcs
        
end



% Scalar error quantities (l1 err)
norm_linf_ref = norm( v_all_CL_c_mp_ref, inf )

err_linf_CL_c_low = norm( v_all_CL_c_low - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;
err_linf_CL_c_med = norm( v_all_CL_c_med - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;
err_linf_CL_c_high = norm( v_all_CL_c_high - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;

err_linf_SH_c_low = norm( v_all_SH_c_low - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;
err_linf_SH_c_high = norm( v_all_SH_c_high - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;

err_linf_PF_c_low = norm( v_all_PF_c_low - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;
err_linf_PF_c_med = norm( v_all_PF_c_med - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;
err_linf_PF_c_high = norm( v_all_PF_c_high - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;

err_linf_ADPTV_c_low = norm( v_all_ADPTV_c_low - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;
err_linf_ADPTV_c_med = norm( v_all_ADPTV_c_med - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;
err_linf_ADPTV_c_high = norm( v_all_ADPTV_c_high - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;

err_linf_PFmp_c_low = norm( v_all_PFmp_c_low - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;
err_linf_PFmp_c_med = norm( v_all_PFmp_c_med - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;
err_linf_PFmp_c_high = norm( v_all_PFmp_c_high - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;

err_linf_DIM_low = norm( v_all_DIM_low - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;
err_linf_DIM_med = norm( v_all_DIM_med - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;
err_linf_DIM_high = norm( v_all_DIM_high - v_all_CL_c_mp_ref, inf ) / norm_linf_ref;




fprintf( 'CLlow: %0.2e\n', err_linf_CL_c_low );
fprintf( 'CLmed: %0.2e\n', err_linf_CL_c_med );
fprintf( 'CLhigh: %0.2e\n', err_linf_CL_c_high );

fprintf( 'SHlow: %0.2e\n', err_linf_SH_c_low );
fprintf( 'SHhigh: %0.2e\n', err_linf_SH_c_high );

fprintf( 'PFlow: %0.2e\n', err_linf_PF_c_low );
fprintf( 'PFmed: %0.2e\n', err_linf_PF_c_med );
fprintf( 'PFhigh: %0.2e\n', err_linf_PF_c_high );

fprintf( 'ADPTVlow: %0.2e\n', err_linf_ADPTV_c_low );
fprintf( 'ADPTVmed: %0.2e\n', err_linf_ADPTV_c_med );
fprintf( 'ADPTVhigh: %0.2e\n', err_linf_ADPTV_c_high );

fprintf( 'PFmp32: %0.2e\n', err_linf_PFmp_c_low );
fprintf( 'PFmp48: %0.2e\n', err_linf_PFmp_c_med );
fprintf( 'PFmp64: %0.2e\n', err_linf_PFmp_c_high );

fprintf( 'DIMlow: %0.2e\n', err_linf_DIM_low );
fprintf( 'DIMmed: %0.2e\n', err_linf_DIM_med );
fprintf( 'DIMhigh: %0.2e\n', err_linf_DIM_high );




figure(1);
plt_CL_low = loglog( v_Nk, v_t_CL_c_low, ':', 'Color', v_col_blue, 'LineWidth', 4.0 );
hold on;
plt_CL_med = loglog( v_Nk, v_t_CL_c_med, '--', 'Color', v_col_blue, 'LineWidth', 4.0 );
plt_CL_high = loglog( v_Nk, v_t_CL_c_high, '-', 'Color', v_col_blue, 'LineWidth', 4.0 );

% plt_PF_low = loglog( v_Nk, v_t_PF_c_low, 'g:', 'LineWidth', 4.0 );
% plt_PF_med = loglog( v_Nk, v_t_PF_c_med, 'g--', 'LineWidth', 4.0 );
% plt_PF_high = loglog( v_Nk, v_t_PF_c_high, 'g-', 'LineWidth', 4.0 );

plt_ADPTV_low = loglog( v_Nk, v_t_ADPTV_c_low, ':', 'Color', v_col_green, 'LineWidth', 4.0 );
plt_ADPTV_med = loglog( v_Nk, v_t_ADPTV_c_med, '--', 'Color', v_col_green, 'LineWidth', 4.0 );
plt_ADPTV_high = loglog( v_Nk, v_t_ADPTV_c_high, '-', 'Color', v_col_green, 'LineWidth', 4.0 );

%plt_PFmp_low = loglog( v_Nk, v_t_PFmp_c_low, 'm-', 'LineWidth', 4.0 );
%plt_PFmp_med = loglog( v_Nk, v_t_PFmp_c_med, 'm-', 'LineWidth', 4.0 );
%plt_PFmp_high = loglog( v_Nk, v_t_PFmp_c_high, 'm-', 'LineWidth', 4.0 );

plt_SH_low = loglog( v_Nk, v_t_SH_c_low, ':', 'Color', v_col_yellow, 'LineWidth', 4.0 );
plt_SH_high = loglog( v_Nk, v_t_SH_c_high, '-', 'Color', v_col_yellow, 'LineWidth', 4.0 );

plt_DIM_low = loglog( v_Nk, v_t_DIM_low, ':', 'Color', v_col_red, 'LineWidth', 4.0 );
plt_DIM_med = loglog( v_Nk, v_t_DIM_med, '--', 'Color', v_col_red, 'LineWidth', 4.0 );
plt_DIM_high = loglog( v_Nk, v_t_DIM_high, '-', 'Color', v_col_red, 'LineWidth', 4.0 );


% % Text labels
% text( 1.5e4, 0.3, 'PF high', 'Fontsize', 18 );



% Ok, the shooting is misnamed in our code, it should be med instead of
% low.
lgnd = legend( [ plt_CL_low, plt_CL_med, plt_CL_high, ... 
                plt_ADPTV_low, plt_ADPTV_med, plt_ADPTV_high, ...
                plt_DIM_low, plt_DIM_med, plt_DIM_high, ...
                plt_SH_low, plt_SH_high ], ...
    sprintf( 'CL low' ), ...
	sprintf( 'CL med' ), ...
    sprintf( 'CL high' ), ...
	sprintf( 'PF low' ), ...
	sprintf( 'PF med' ), ...
	sprintf( 'PF high' ), ...
	sprintf( 'DIM low' ), ...
	sprintf( 'DIM med' ), ...
	sprintf( 'DIM high' ), ... 
	sprintf( 'SH med' ), ...
	sprintf( 'SH high' ) );


% lgnd = legend( [ plt_CL_low, plt_CL_med, plt_CL_high, ... 
%                 plt_ADPTV_low, plt_ADPTV_med, plt_ADPTV_high, ...
%                 plt_SH_low, plt_SH_high, ...
%                 plt_DIM_low, plt_DIM_med, plt_DIM_high ], ...
%     sprintf( 'CL low; err=%0.0e.', err_linf_CL_c_low ), ...
% 	sprintf( 'CL med; err=%0.0e.', err_linf_CL_c_med ), ...
%     sprintf( 'CL high; err=%0.0e.', err_linf_CL_c_high ), ...
% 	sprintf( 'PF low; err=%0.0e.', err_linf_ADPTV_c_low ), ...
% 	sprintf( 'PF med; err=%0.0e.', err_linf_ADPTV_c_med ), ...
% 	sprintf( 'PF high; err=%0.0e.', err_linf_ADPTV_c_high ), ...
% 	sprintf( 'SH med; err=%0.0e.', err_linf_SH_c_low ), ...
% 	sprintf( 'SH high; err=%0.0e.', err_linf_SH_c_high ), ...
% 	sprintf( 'DIM low; err=%0.0e.', err_linf_DIM_low ), ...
% 	sprintf( 'DIM med; err=%0.0e.', err_linf_DIM_med ), ...
% 	sprintf( 'DIM high; err=%0.0e.', err_linf_DIM_high ) );



% lgnd = legend( [ plt_CL_64, plt_DIM_high, plt_PF_high ], ...
%     sprintf( 'Collocation', err_linf_CL_c_high ),     ...
% 	sprintf( 'DIM', err_linf_DIM_high ), ...
% 	sprintf( 'Path-following', err_linf_PF_c_high ) );




xlim( [ 0 v_Nk(end) ] );
xticks( [ 1e2 1e3 1e4 1e5 1e6 1e7 ] );
ylim( [ 1e-3 1e3 ] );
yticks( [ 1e-3 1e-2 1e-1 1e-0 1e1 1e2 1e3 ] );
% ylim( [ 1e-3 1e7 ] );
% yticks( [ 1e-3 1e-2 1e-1 1e-0 1e1 1e2 1e3 1e4 1e5 1e6 1e7 ] );


set(lgnd,'Location','southoutside');
set(lgnd,'Interpreter','latex');
set(lgnd,'FontSize', 18 );
set(lgnd,'NumColumns', 4 );

set(gca,'fontsize', 22 );
set(gca,'TickLabelInterpreter', 'latex');

xlabel( 'Number of $k$ query points', 'Interpreter', 'Latex', 'FontSize', 22 );
ylabel( 'Time (s)', 'Interpreter', 'Latex', 'FontSize', 22 );

grid on;

hold off;
    

set(gcf,'PaperSize', [ 24 26 ] );
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\perf_std', '-dpdf' );


% set(gcf,'PaperSize', [ 20 32 ] );
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'PaperPositionMode','auto');
% set(gcf,'PaperUnits','normalized');
% set(gcf,'PaperPosition', [0 0 1 1]);
% print( 'figures\perf_std', '-dpdf' );

    
end
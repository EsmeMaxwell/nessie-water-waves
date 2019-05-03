function fn_ww__analyse_re__acc_tests_c__red(  )
%fn_ww__analyse_re__acc_tests_c__red: Plot eig accuracy and berr+cond results for Rayleigh
% 
%   fn_ww__analyse_re__acc_tests_c__red() 
% 
% Uses precalculated data on eigenvalue accuracy, backward error, and
% condition numbers to generate two figures.
% 
% TAGS: SISCPFLIB
%
% REQUIRES
% 
% Files: st_acc_tests__ref and st_acc_tests__data% 
% 
% See also
%   fn_ww__sim_re__acc_tests_c__red()


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



s_datapath = 'sisc_data';


% Load reference set
load( sprintf( '%s\\st_acc_tests__ref', s_datapath ) );

% Load dataset
load( sprintf( '%s\\st_acc_tests__data', s_datapath ) );


norm_ref = norm( st_ref.v_c_CL_ref, inf );

v_c_err_CL = zeros( st_data.NNz, 1 );
v_c_err_PF = zeros( st_data.NNz, 1 );
v_c_err_PFmp = zeros( st_data.NNz, 1 );
v_c_err_DIM = zeros( st_data.NNz, 1 );

v_berr_CL_inf = zeros( st_data.NNz, 1 );   % inf norm of v_berr
v_cond_CL_inf = zeros( st_data.NNz, 1 );   % inf norm of v_cond
v_berr_PF_inf = zeros( st_data.NNz, 1 );   % inf norm of v_berr
v_cond_PF_inf = zeros( st_data.NNz, 1 );   % inf norm of v_cond


for lp_Nz=1:st_data.NNz
    
    % Calc errors (normwise, l2)
    v_c_err_CL( lp_Nz ) = norm( st_data.ca_v_c_CL{lp_Nz,1} - st_ref.v_c_CL_ref, inf ) / norm_ref;
    v_c_err_PF( lp_Nz ) = norm( st_data.ca_v_c_PF{lp_Nz,1} - st_ref.v_c_CL_ref, inf ) / norm_ref;
    v_c_err_PFmp( lp_Nz ) = norm( st_data.ca_v_c_PFmp{lp_Nz,1} - st_ref.v_c_CL_ref, inf ) / norm_ref;
    v_c_err_DIM( lp_Nz ) = norm( st_data.ca_v_c_DIM{lp_Nz,1} - st_ref.v_c_CL_ref, inf ) / norm_ref;
    
    % Collcate berr and cond
    v_berr_CL_inf( lp_Nz ) = norm( st_data.ca_st_err_CL{lp_Nz}.v_berr, inf );
    v_cond_CL_inf( lp_Nz ) = norm( st_data.ca_st_err_CL{lp_Nz}.v_cond, inf );
    v_berr_PF_inf( lp_Nz ) = norm( st_data.ca_st_err_PF{lp_Nz}.v_berr, inf );
    v_cond_PF_inf( lp_Nz ) = norm( st_data.ca_st_err_PF{lp_Nz}.v_cond, inf );
    
end


% Temporary fix for bad line style, remove every 2nd data entry
st_data.v_Nz = st_data.v_Nz(1:2:end);
v_c_err_CL = v_c_err_CL(1:2:end);
v_c_err_PF = v_c_err_PF(1:2:end);
v_c_err_PFmp = v_c_err_PFmp(1:2:end);
v_c_err_DIM = v_c_err_DIM(1:2:end);

v_berr_CL_inf = v_berr_CL_inf(1:2:end);
v_cond_CL_inf = v_cond_CL_inf(1:2:end);
v_berr_PF_inf = v_berr_PF_inf(1:2:end);
v_cond_PF_inf = v_cond_PF_inf(1:2:end);


% 
% 
% figure(1);
% loglog( st_data.v_Nz, v_c_err_CL, '-', 'Color', v_col_blue, 'LineWidth', 5.0 );
% hold on;
% loglog( st_data.v_Nz, v_c_err_PF, '-.', 'Color', v_col_green, 'LineWidth', 5.0 ); 
% loglog( st_data.v_Nz, v_c_err_PFmp, '--', 'Color', v_col_purple, 'LineWidth', 5.0 );
% loglog( st_data.v_Nz, v_c_err_DIM, ':', 'Color', v_col_red, 'LineWidth', 5.0 );
% hold off;
% 
% text( 70, 1e-14, 'PFmp', 'Fontsize', 28, 'Interpreter', 'Latex' );
% text( 110, 1e-10, 'PF', 'Fontsize', 28, 'Interpreter', 'Latex' );
% text( 75, 1e-7, 'CL', 'Fontsize', 28, 'Interpreter', 'Latex' );
% text( 120, 1e-4, 'DIM', 'Fontsize', 28, 'Interpreter', 'Latex' );
% 
% 
% 
% lgnd = legend( 'CL', 'PF', 'PFmp', 'DIM' );
% 
% set(lgnd,'Location','southwest');
% set(lgnd,'Interpreter','latex');
% set(lgnd,'FontSize', 32 );
% 
% 
% set(gca,'fontsize', 32 );
% set(gca,'TickLabelInterpreter', 'latex');
% 
% xlabel( '$N_z$', 'Interpreter', 'Latex', 'FontSize', 32 );
% ylabel( 'Error $\epsilon$', 'Interpreter', 'Latex', 'FontSize', 32 );
% %ylabel( 'Error $\epsilon$', 'Interpreter', 'Latex', 'FontSize', 28 );
% 
% xlim( [ 20 160 ] );
% ylim( [ 1e-16 1e-1 ] );
% xtickangle(90);
% xticks( 20:20:160 );
% yticks( [ 1e-16 1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 ] );
% 
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'PaperPositionMode','auto');
% set(gcf,'PaperUnits','normalized');
% set(gcf,'PaperPosition', [0 0 1 1]);
% print( 'figures\eig_accuracy', '-dpdf' );
% 










figure(2);
loglog( st_data.v_Nz, v_berr_CL_inf, '--', 'Color', v_col_blue, 'LineWidth', 5.0 );
hold on;
loglog( st_data.v_Nz, v_cond_CL_inf, ':', 'Color', v_col_blue, 'LineWidth', 5.0 );
loglog( st_data.v_Nz, v_cond_CL_inf .* v_berr_CL_inf, '-', 'Color', v_col_blue, 'LineWidth', 5.0 );

loglog( st_data.v_Nz, v_berr_PF_inf, '--', 'Color', v_col_purple, 'LineWidth', 5.0 );
loglog( st_data.v_Nz, v_cond_PF_inf, ':', 'Color', v_col_purple, 'LineWidth', 5.0 );
loglog( st_data.v_Nz, v_cond_PF_inf .* v_berr_PF_inf, '-', 'Color', v_col_purple, 'LineWidth', 5.0 );


text( 30, 1e8, 'PF cond', 'Fontsize', 28, 'Interpreter', 'Latex' );
text( 20, 3.5e2, 'CL cond', 'Fontsize', 28, 'Interpreter', 'Latex' );

text( 70, 1e-5, 'CL berr $\times$ cond', 'Fontsize', 28, 'Interpreter', 'Latex' );
text( 70, 2.5e-11, 'PG berr $\times$ cond', 'Fontsize', 28, 'Interpreter', 'Latex' );

text( 35, 2.5e-13, 'CL berr', 'Fontsize', 28, 'Interpreter', 'Latex' );
text( 35, 2.5e-18, 'PF berr', 'Fontsize', 28, 'Interpreter', 'Latex' );


% lgnd = legend( 'CL berr', 'CL cond', 'CL berr $\times$ cond', 'PF berr', 'PF cond', 'PF berr $\times$ cond' );
% 
% set(lgnd,'Interpreter','latex');
% set(lgnd,'FontSize', 26 );
% set(lgnd,'NumColumns', 2 );
% set(lgnd,'Position', [ 0.1 0.5 1.0469 0.3090] );


set(gca,'fontsize', 28 );
set(gca,'TickLabelInterpreter', 'latex');

xlabel( '$N_z$', 'Interpreter', 'Latex', 'FontSize', 28 );
ylabel( '$l^{\infty}$ of qty over $k$', 'Interpreter', 'Latex', 'FontSize', 28 );


xlim( [ 0 160 ] );
ylim( [ 1e-20 1e+12 ] );
xtickangle(90);
xticks( 0:20:160 );
yticks( [ 1e-20 1e-16 1e-12 1e-8 1e-4 1e0 1e4 1e8 1e12 ] );


hold off;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print( 'figures\eig_berr_cond', '-dpdf' );




end
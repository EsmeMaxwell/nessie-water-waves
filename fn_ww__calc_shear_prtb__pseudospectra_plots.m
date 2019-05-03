function fn_ww__calc_shear_prtb__pseudospectra_plots( Ngrid )
%fn_ww__calc_shear_prtb__pseudospectra_plots: Calc pseudospectra for the err-in-shear-profile problem
%
%   fn_ww__calc_shear_prtb__pseudospectra_plots( Ngrid ) 
% 
% 
% TAGS: WWERRINSHEAR
%



set(0,'defaulttextinterpreter','latex');



% EXP PLOT-----------------------------------------------------------------

phy_h = 10;
h = phy_h;

Nz = 36;
k=pi;

% Setup parameter space
[ st_p ] = fn_ww__setup__param_std__re_cl(  );

% Setup profile and update parameter set    
[ st_fn_shear_exp, st_p ] = fn_ww__setup__shear_fn__nondim_exp( st_p, 0.05, phy_h, 1 );
[ st_r_shear_exp ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear_exp, st_p );

% Setup differentiation matrix, etc.
[ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( Nz, 1 );
[ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p );

% Get shear data
[ st_v_shear ] = fn_ww__util__get_r_shear_data( st_Dn, st_r_shear_exp, st_p );

% Calculate reference v_c
st_v_shear = struct;
st_v_shear.v_U = st_fn_shear_exp.fn_U( st_Dn.v_zm );
st_v_shear.v_dU = st_fn_shear_exp.fn_dU( st_Dn.v_zm );
st_v_shear.v_ddU = st_fn_shear_exp.fn_ddU( st_Dn.v_zm );





% %POWERLAW --------------------------------------------------------
% 
% Fr2 = 0.05;
% phy_h = 10;
% h = phy_h;
% 
% Nz = 48;
% k=pi;
% 
% % Setup parameter space
% [ st_p ] = fn_ww__setup__param_std__re_cl(  );
% 
% % Setup profile and update parameter set
% phy_U0 = sqrt( Fr2 * st_p.phy_g * phy_h );
% [ st_fn_shear_pwr, st_p ] = fn_ww__setup__shear_fn__nondim_powerlaw( st_p, phy_U0, phy_h );
% [ st_r_shear_pwr ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear_pwr, st_p );
% 
% % Setup differentiation matrix, etc.
% [ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( Nz, 1 );
% [ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p );
% 
% % Get shear data
% [ st_v_shear ] = fn_ww__util__get_r_shear_data( st_Dn, st_r_shear_pwr, st_p );    
% 
% % % Calculate reference v_c
% % st_v_shear = struct;
% % st_v_shear.v_U = st_fn_shear_pwr.fn_U( st_Dn.v_zm );
% % st_v_shear.v_dU = st_fn_shear_pwr.fn_dU( st_Dn.v_zm );
% % st_v_shear.v_ddU = st_fn_shear_pwr.fn_ddU( st_Dn.v_zm );

    
    
    
% 
% 
% % COLOMBIA -----------------------------------------
% 
% phy_h = 20;
% 
% h = phy_h;
% 
% Nz = 48;
% k=pi;
% 
% % Setup parameter space
% [ st_p ] = fn_ww__setup__param_std__re_cl(  );
% 
% % Setup profile and update parameter set    
% [ st_fn_shear_poly, st_p ] = fn_ww__setup__shear_fn__nondim_columbia_poly( st_p );
% [ st_r_shear_poly ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear_poly, st_p );
% 
% % Setup differentiation matrix, etc.
% [ st_Dn ] = fn_ww__setup__diffmtrx__WR_poldif( Nz, 1 );
% [ st_Dn ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p );
% 
% % Get shear data
% [ st_v_shear ] = fn_ww__util__get_r_shear_data( st_Dn, st_r_shear_poly, st_p );
% 
% % Calculate reference v_c
% st_v_shear = struct;
% st_v_shear.v_U = st_fn_shear_poly.fn_U( st_Dn.v_zm );
% st_v_shear.v_dU = st_fn_shear_poly.fn_dU( st_Dn.v_zm );
% st_v_shear.v_ddU = st_fn_shear_poly.fn_ddU( st_Dn.v_zm );
% 
% 

    
    
    
    






% U matrices
a_U = diag( st_v_shear.v_U );
a_dU = diag( st_v_shear.v_dU );
a_ddU = diag( st_v_shear.v_ddU );
a_I = eye( size( a_U ) );

% Precomputations
a_A0_precomp = (  a_U * st_Dn.a_D2m - a_ddU  );
a_A2 = 0 * a_I;

% FS
a_A0_FS = ( a_U^2 * st_Dn.a_Dm - a_dU * a_U - ( 1 / st_p.Fr2 ) * a_I );
a_A1_FS = ( -2 * a_U * st_Dn.a_Dm + a_dU );
a_A2_FS = st_Dn.a_Dm;


    
% Main equations
a_A0 = ( a_A0_precomp - k^2 * a_U );
a_A1 = -st_Dn.a_D2m + k^2 * a_I;

a_A0(1,:) = a_A0_FS(1,:);
a_A1(1,:) = a_A1_FS(1,:);
a_A2(1,:) = a_A2_FS(1,:);

a_A0(end,:) = 0 * a_A0(end,:);
a_A1(end,:) = 0 * a_A1(end,:);
a_A2(end,:) = 0 * a_A2(end,:);




% Calc eigenvalues
[ v_eigs ] = polyeig( a_A0, a_A1, a_A2 );
%v_eigs = diag( a_eig );
v_eigs( ~isfinite(v_eigs) ) = [];
v_eigs_sorted = sort( v_eigs, 'descend', 'ComparisonMethod','real' );

% Calc scalings
% alpha_2 = norm( a_A2, 2 )
% alpha_1 = norm( a_A1, 2 )
% alpha_0 = norm( a_A0, 2 )
alpha_2 = 1;
alpha_1 = 1;
alpha_0 = 1;

% % EXP
% min_re = -8;
% max_re = 12;
% min_im = -10;
% max_im = 10;

min_re = -10;
max_re = 10;
min_im = -10;
max_im = 10;

nx = Ngrid;
ny = Ngrid;
v_x = linspace( min_re, max_re, nx );
v_y = linspace( min_im, max_im, ny );

Z = zeros( nx, ny );
for lp_x=1:nx
    lp_x
    for lp_y=1:ny
        
        % Relevant c
        c = v_x(lp_x) + 1i*v_y(lp_y);

        % Evaluate p poly
        p = alpha_2 * abs(c)^2 + alpha_1 * abs(c) + alpha_0;
        
        % Do resolvent calc
        a_P = ( c^2 * a_A2 + c * a_A1 + a_A0 );
        a_P = a_P(1:(end-1),1:(end-1));
        a_Pinv = inv( a_P );        
        Z(lp_y,lp_x) = norm( a_Pinv, 2 );
        
        % Finalise
        Z(lp_y,lp_x) = Z(lp_y,lp_x) * p;
               
    end
end



% st_PS_c = struct;
% st_PS_c.v_eigs_sorted = v_eigs_sorted;
% st_PS_c.nx = nx;
% st_PS_c.ny = ny;
% st_PS_c.grid_ext = grid_ext;
% st_PS_c.v_x = v_x;
% st_PS_c.v_y = v_y;
% st_PS_c.Z = Z;
% save( 'st_PS_c', 'st_PS_c' );

% Fontsizes
fontsz_labels = 26;
fontsz_axes = 26;
fontsz_title = 26;




%v_levels = 0.05:0.05:3;  %EXP
v_levels = 0.1:0.1:2.5;  % Powerlaw

h_fig = figure(1);
colormap bone;
[ a_M, ob_c ] = contour( v_x, v_y, 1./Z, v_levels, 'ShowText', 'on', 'LabelSpacing', 288 );
%[ a_M, ob_c ] = contour( v_x, v_y, Z, 40, 'ShowText', 'on', 'LabelSpacing', 288 );
ob_c.LineWidth = 1.5;
hold on;
plot(real(v_eigs_sorted(2:end)),imag(v_eigs_sorted(2:end)), 'b*', 'MarkerSize', 16, 'LineWidth', 3.0 );
plot(real(v_eigs_sorted(1)),imag(v_eigs_sorted(1)), 'r*', 'MarkerSize', 24, 'LineWidth', 4.0 );
xlim( [ min_re max_re ] );
ylim( [ min_im max_im ] );
set(gca,'fontsize', fontsz_axes );
set(gca,'TickLabelInterpreter', 'latex');        
xlabel( '$\Re\, \tilde{c}$', 'Interpreter', 'latex', 'FontSize', fontsz_labels );
ylabel( '$\Im\, \tilde{c}$', 'Interpreter', 'latex', 'FontSize', fontsz_labels );
hold off;


set(gcf,'PaperOrientation','portrait');  %Changed
set(gcf,'PaperSize', [ 20 20 ] );      
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);

set( h_fig, 'Color', 'none' );
set( h_fig, 'InvertHardCopy', 'Off' );        
print( 'figures_shear_prtb\fig__ps_exp_kPi_Fr0_05', '-dpdf' );
%print( 'figures_shear_prtb\fig__ps_powerlaw_kPi_Fr0_05', '-dpdf' );
%print( 'figures_shear_prtb\fig__ps_columbia_kPi_Fr0_05', '-dpdf' );





% figure(2);
% [ a_x, a_y ] = meshgrid( v_x, v_y );
% surface( a_x, a_y, Z );


% min( min( Z ) )
% max( max( Z ) ) 
% 
% v_eigs_sorted

end
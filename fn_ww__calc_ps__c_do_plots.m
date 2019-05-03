function fn_ww__calc_ps__c_do_plots( st_Dn, k, st_r_shear, Ngrid, st_p )
%fn_ww__calc_ps__c_do_plots: Calc & plot pseudospectrum for given shear profile
% 
%   fn_ww__calc_ps__c_do_plots( st_Dn, k, st_v_shear, Ngrid, st_p )
% 
% Calculates and plots the pseudospectrum for given shear profile
% (polynomial pseudospectrum)
% 


% TODO change to return whole eval set and eigenvectors

% TODO check whether error is due to mix-up with a,b interval

N = numel( st_Dn.v_z0 );


% Get shear data
[ st_v_shear ] = fn_ww__util__get_r_shear_data( st_Dn, st_r_shear, st_p );    

% U matrices
a_U = diag( st_v_shear.v_U );
a_dU = diag( st_v_shear.v_dU );
a_ddU = diag( st_v_shear.v_ddU );
a_I = eye( size( a_U ) );

% Precalc
[ a_A0_precomp, a_A2, a_A0_FS, a_A1_FS, a_A2_FS ] = fn_ww__calc_re__cl__mtrxs__precalc_c( st_Dn, a_U, a_dU, a_ddU, a_I, st_p );


% Core equations
[ a_A0, a_A1, a_A2 ] = fn_ww__calc_re__cl__mtrxs__core_c( st_Dn, a_U, a_A0_precomp, a_A2, a_I, a_A0_FS, a_A1_FS, a_A2_FS, k, st_p );

% Calc eigenvalues
[ v_eigs ] = polyeig( a_A0, a_A1, a_A2 );
%v_eigs = diag( a_eig );
v_eigs( ~isfinite(v_eigs) ) = [];
v_eigs_sorted = sort( v_eigs, 'descend', 'ComparisonMethod','real' );

% Scalings
alpha_2 = 1;
alpha_1 = 1;
alpha_0 = 1;


% Get min, max
% min_re = min( real( v_eigs_sorted(1:12) ) )
% max_re = max( real( v_eigs_sorted(1:12) ) )
% min_im = min( imag( v_eigs_sorted(1:12) ) )
% max_im = max( imag( v_eigs_sorted(1:12) ) )

min_re = -12;
max_re = 20;
min_im = -16;
max_im = 16;


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

figure(1);
contour( v_x, v_y, Z, 40 );
hold on;
plot(real(v_eigs_sorted(2:end)),imag(v_eigs_sorted(2:end)), 'b*', 'MarkerSize', 12, 'LineWidth', 2.0 );
plot(real(v_eigs_sorted(1)),imag(v_eigs_sorted(1)), 'r*', 'MarkerSize', 12, 'LineWidth', 2.0 );
xlim( [ min_re max_re ] );
ylim( [ min_im max_im ] );
xlabel( 'Re', 'Interpreter', 'latex', 'FontSize', 16 );
ylabel( 'Im', 'Interpreter', 'latex', 'FontSize', 16 );
hold off;

figure(2);
[ a_x, a_y ] = meshgrid( v_x, v_y );
surface( a_x, a_y, Z );


min( min( Z ) )
max( max( Z ) ) 

v_eigs_sorted

end
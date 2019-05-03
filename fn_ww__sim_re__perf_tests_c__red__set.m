function  fn_ww__sim_re__perf_tests_c__red__set( )
%fn_ww__sim_re__perf_tests_c__red__set: Simulation (batch mgr) for CL,DIM,PF perf tests
% 
% 	fn_ww__sim_re__perf_tests_c__red__set()
% 
% Performs the performance tests, calling
% fn_ww__sim_re__perf_tests_c__red() several times.
% 
% TAGS: SISCPFLIB
% 
% See also
%   fn_ww__sim_re__perf_tests_c__red(),
%   fn_ww__sim_re__perf_tests_c__red__parfor()



k_min = 0.025;
k_max = 250;



[ st_p ] = fn_ww__setup__param_std__re_cl(  );
Fr2 = 0.05; phy_h = 20;
[ st_fn_shear, st_p ] = fn_ww__setup__shear_fn__nondim_cospwr( st_p, sqrt( Fr2 * st_p.phy_g * phy_h ), phy_h, 2, 4 * pi, 1, 0.5 );



% % Do ref calculations here
% fprintf( '\n\n\n[ v_Nk = 100 ]\n' );
% v_Nk = [ 100 ];
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk100_ref__k0_025_250__cospwr', 1, 1, 1, 1, 0, st_p );
% 
% 
% % No ref calcs, still do CL
% fprintf( '\n\n\n[ v_Nk = 300:50:1000 ]\n' );
% v_Nk = 300:50:1000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk300_50_1000__k0_025_250__cospwr', 0, 1, 1, 1, 1, st_p );
% 
% 
% % No ref calcs, still do CL
% fprintf( '\n\n\n[ v_Nk = 1125:125:2500 ]\n' );
% v_Nk = 1125:125:2500;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk1125_125_2500__k0_025_250__cospwr', 0, 1, 1, 1, 1, st_p );
% 
% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 3000:250:5000 ]\n' );
% v_Nk = 3000:250:5000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk3000_250_5000__k0_025_250__cospwr', 0, 1, 1, 1, 1, st_p );


% % No ref calcs, no CL (memory issues arise here)
% fprintf( '\n\n\n[ v_Nk = 6000:500:7500 ]\n' );
% v_Nk = 6000:500:7500;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk6000_500_7500__k0_025_250__cospwr', 0, 1, 1, 1, 1, st_p );
% 
% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 8000:1000:10000 ]\n' );
% v_Nk = 8000:1000:10000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk8000_1000_10000__k0_025_250__cospwr', 0, 1, 1, 0, 1, st_p );


% No ref calcs (only one repetition)
fprintf( '\n\n\n[ v_Nk = 11000:1000:12000 ]\n' );
v_Nk = 11000:1000:12000;
fn_ww__sim_re__perf_tests_c__red( v_Nk, 1, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk11000_1000_12000__k0_025_250__cospwr', 0, 1, 1, 0, 1, st_p );

% No ref calcs (only one repetition)
fprintf( '\n\n\n[ v_Nk = 13000:1000:14000 ]\n' );
v_Nk = 13000:1000:14000;
fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk13000_1000_14000__k0_025_250__cospwr', 0, 1, 1, 0, 1, st_p );

% No ref calcs (only one repetition)
fprintf( '\n\n\n[ v_Nk = 15000:1000:16000 ]\n' );
v_Nk = 15000:1000:16000;
fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk15000_1000_16000__k0_025_250__cospwr', 0, 1, 1, 0, 1, st_p );

% No ref calcs (only one repetition)
fprintf( '\n\n\n[ v_Nk = 17000:1000:18000 ]\n' );
v_Nk = 17000:1000:18000;
fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk17000_1000_18000__k0_025_250__cospwr', 0, 1, 1, 0, 1, st_p );

% No ref calcs (only one repetition)
fprintf( '\n\n\n[ v_Nk = 19000:1000:20000 ]\n' );
v_Nk = 19000:1000:20000;
fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk19000_1000_20000__k0_025_250__cospwr', 0, 1, 1, 0, 1, st_p );



% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 21000:1000:30000 ]\n' );
% v_Nk = 21000:1000:30000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk21000_1000_30000__k0_025_250__cospwr', 0, 0, 1, 0, 1, st_p );
% 
% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 32500:2500:50000 ]\n' );
% v_Nk = 32500:2500:50000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk32500_2500_50000__k0_025_250__cospwr', 0, 0, 1, 0, 1, st_p );
% 
% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 55000:5000:100000 ]\n' );
% v_Nk = 55000:5000:100000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk55000_5000_100000__k0_025_250__cospwr', 0, 0, 1, 0, 1, st_p );
% 
% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 150000:10000:200000 ]\n' );
% v_Nk = 150000:10000:200000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk150000_10000_200000__k0_025_250__cospwr', 0, 0, 1, 0, 1, st_p );
% 
% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 250000:50000:500000 ]\n' );
% v_Nk = 250000:50000:500000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk250000_50000_500000__k0_025_250__cospwr', 0, 0, 0, 0, 1, st_p );
% 
% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 600000:100000:1000000 ]\n' );
% v_Nk = 600000:100000:1000000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk600000_100000_1000000__k0_025_250__cospwr', 0, 0, 0, 0, 1, st_p );
% 
% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 1100000:100000:1500000 ]\n' );
% v_Nk =1100000:100000:1500000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk1100000_100000_1500000__k0_025_250__cospwr', 0, 0, 0, 0, 1, st_p );
% 
% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 1600000:100000:2000000 ]\n' );
% v_Nk =1600000:100000:2000000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk1600000_100000_2000000__k0_025_250__cospwr', 0, 0, 0, 0, 1, st_p );
% 
% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 2100000:100000:2500000 ]\n' );
% v_Nk = 2100000:100000:2500000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk2100000_100000_2500000__k0_025_250__cospwr', 0, 0, 0, 0, 1, st_p );
% 
% 
% % No ref calcs, no CL
% fprintf( '\n\n\n[ v_Nk = 2600000:100000:3000000 ]\n' );
% v_Nk = 2600000:100000:3000000;
% fn_ww__sim_re__perf_tests_c__red( v_Nk, 5, st_fn_shear, k_min, k_max, 'sim_perf_c_set29__Nk2600000_100000_3000000__k0_025_250__cospwr', 0, 0, 0, 0, 1, st_p );
% 
% 

end
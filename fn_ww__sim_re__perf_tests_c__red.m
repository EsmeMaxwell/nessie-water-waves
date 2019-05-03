function  fn_ww__sim_re__perf_tests_c__red( v_Nk, Nrep, st_fn_shear, k_min, k_max, s_filename, b_do_ref, b_do_CL, b_do_DIMmed, b_do_DIMhigh, b_no_save, st_p )
%fn_ww__sim_re__perf_tests_c__red: Simulation for CL,DIM,PF perf tests
%
%   fn_ww__sim_re__perf_tests_c__red( v_Nk_Ntheta, Nrep, st_CL_fn_shear_x, st_CL_fn_shear_y, k_min, k_max, s_filename, b_do_CL, b_no_save, st_p )
% 
% Runs the performance tests for for algorithms in the 1d context. Each
% test is performed Nrep times and the lowest result is saved; this should
% help ensure that caches are warm before starting each test.
%
% This is usually called in a sequential manner from
% fn_ww__sim_re__perf_tests_c__red__set(), which ensures output files do
% not exceed MATLAB's 2Gb limit or memory limits.
%
% TAGS: SISCPFLIB
% 
% INPUT
% 
% v_Nk : vector with required number of k elements for each test
% 
% Nrep : The number of repetitions for each test.
% 
% st_CL_fn_shear : The shear profile to use for the tests
% 
% k_min, k_max : k interval to use
% 
% s_filename : Output data filename
% 
% b_do_CL : Flag for whether to perform the collocation tests. Anything
% nonzero is true.
% 
% b_no_save : Flag to determine if we save the actual test point data (we
% always save the time data). If nonzero then we do not, this saves an
% immense amount of space but an error estimate cannot be calculated using
% these tests.
% 
% 
% See also
%   fn_ww__sim_re__perf_tests_c__red__parfor(),
%   fn_ww__sim_re__perf_tests_c__red__set()


NNk = numel( v_Nk );


% Setup for CL
[ st_Dn_mp_ref ] = fn_ww__setup__diffmtrx_mp__WR_chebdif( 256 );
[ st_Dn_mp_ref ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_mp_ref, st_p );

[ st_Dn_std_low ] = fn_ww__setup__diffmtrx__WR_poldif( 32, 1 );
[ st_Dn_std_med ] = fn_ww__setup__diffmtrx__WR_poldif( 52, 1 );
[ st_Dn_std_high ] = fn_ww__setup__diffmtrx__WR_poldif( 64, 1 );
[ st_Dn_std_low ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_std_low, st_p );
[ st_Dn_std_med ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_std_med, st_p );
[ st_Dn_std_high ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_std_high, st_p );

[ st_Dn_adptv_low ] = fn_ww__setup__diffmtrx__WR_poldif( 22, 1 );
[ st_Dn_adptv_med ] = fn_ww__setup__diffmtrx__WR_poldif( 32, 1 );
[ st_Dn_adptv_high ] = fn_ww__setup__diffmtrx__WR_poldif( 44, 1 );
[ st_Dn_adptv_low ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_adptv_low, st_p );
[ st_Dn_adptv_med ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_adptv_med, st_p );
[ st_Dn_adptv_high ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_adptv_high, st_p );

% [ st_Dn_adptv_dst ] = fn_ww__setup__diffmtrx__WR_poldif( 160, 1 );
% [ st_Dn_adptv_dst ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_adptv_dst, st_p );

[ st_Dn_mp_low ] = fn_ww__setup__diffmtrx_mp__WR_chebdif( 32 );
[ st_Dn_mp_med ] = fn_ww__setup__diffmtrx_mp__WR_chebdif( 52 );
[ st_Dn_mp_high ] = fn_ww__setup__diffmtrx_mp__WR_chebdif( 64 );
[ st_Dn_mp_low ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_mp_low, st_p );
[ st_Dn_mp_med ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_mp_med, st_p );
[ st_Dn_mp_high ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_mp_high, st_p );

[ st_r_shear ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear, st_p );


% Setup for DIM
st_fn_shear_DIM = struct;
st_fn_shear_DIM.fn_Ux = st_fn_shear.fn_phy_U;
st_fn_shear_DIM.fn_dUx = st_fn_shear.fn_phy_dU;
st_fn_shear_DIM.fn_ddUx = st_fn_shear.fn_phy_ddU;
st_fn_shear_DIM.fn_Uy = @(z) 0*z;
st_fn_shear_DIM.fn_dUy = @(z) 0*z;
st_fn_shear_DIM.fn_ddUy = @(z) 0*z;



% Package up necessary data to make fn non-nested
st_data = struct;
st_data.k_min = k_min;
st_data.k_max = k_max;

st_data.st_Dn_mp_ref = st_Dn_mp_ref;

st_data.st_Dn_std_low = st_Dn_std_low;
st_data.st_Dn_std_med = st_Dn_std_med;
st_data.st_Dn_std_high = st_Dn_std_high;

st_data.st_Dn_adptv_low = st_Dn_adptv_low;
st_data.st_Dn_adptv_med = st_Dn_adptv_med;
st_data.st_Dn_adptv_high = st_Dn_adptv_high;

%st_data.st_Dn_adptv_dst = st_Dn_adptv_dst;

st_data.st_Dn_mp_low = st_Dn_mp_low;
st_data.st_Dn_mp_med = st_Dn_mp_med;
st_data.st_Dn_mp_high = st_Dn_mp_high;


st_data.st_r_shear = st_r_shear;
st_data.st_fn_shear_DIM = st_fn_shear_DIM;




% Prepare (moving v_k generation here away from the parfor fn). This means
% we can precompute the reference mp runs.
ca_v_k_q = cell( 1, NNk );
for lp_Nk=1:NNk

    % Use quadratically distributed
    ca_v_k_q{lp_Nk} = fn_ww__util__create_k_vec( st_data.k_min, st_data.k_max, v_Nk(lp_Nk), 3, 0 );
    
end



% If we need reference data, check for saved and if not calc.
if ( b_do_ref )
    
    % Create filename
    s_ref_filename = sprintf( '%s__refdata.mat', s_filename );

    if ( isfile( s_ref_filename ) )

        % Load save ref data
        fprintf( 'Loading ref data...' );
        load( s_ref_filename );
        fprintf( '... done.\n' );       
        
    else        
        
        error( 'Should not be doing this' );
        
        % Calculate reference data (no need for parfor here as Advanpix
        % does multithread)
        fprintf( 'Calculating ref data...\n' );
        ca_st_ref_results = cell( NNk, 1 );
        for lp_Nk=1:NNk
        
            fprintf( '---[ lp_Nk=%d ]\n', lp_Nk );
            [ st_p_ref ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_mp', true, 'bp_disp_update', true ) );           
            [ v_c_CL_mp_ref, ~ ] = fn_ww__calc_re__cl__red_c( st_Dn_mp_ref, ca_v_k_q{lp_Nk}, st_r_shear, st_p_ref );
            ca_st_ref_results{lp_Nk}.v_c_CL_mp_ref = v_c_CL_mp_ref;
            
        end
        save( s_ref_filename, 'ca_st_ref_results' );    
        fprintf( '... done.\n' );        
        
    end
    
else
    
    % We don't need the calculated data for ref but we do need to create
    % empty dataset
    ca_st_ref_results = cell( NNk, 1 );
    for lp_Nk=1:NNk
        ca_st_ref_results{lp_Nk}.v_c_CL_mp_ref = -1;
    end
    
end






parpool(3);
% Parfor loop so have to do annoying stuff to make sure it can index
ca_st_results = cell( NNk, 1 );
parfor lp_Nk=1:NNk
% warning( 'NOT USING PARFOR!' );
% for lp_Nk=1:NNk
    fprintf( '----[ Starting %d of %d ]------\n', lp_Nk, NNk );
    
    % Do the calculation(s)
    [ ca_st_results{lp_Nk} ] = fn_ww__sim_re__perf_tests_c__red__parfor( st_data, ca_v_k_q{lp_Nk}, Nrep, b_do_ref, b_do_CL, b_do_DIMmed, b_do_DIMhigh, b_no_save, st_p );  
        
    fprintf( '----[ End %d of %d ]------\n', lp_Nk, NNk );    
end
delete(gcp);




% Reassemble
v_t_CL_c_low = zeros( NNk, 1 );
v_t_CL_c_med = zeros( NNk, 1 );
v_t_CL_c_high = zeros( NNk, 1 );

v_c_SH_low = zeros( NNk, 1 );
v_c_SH_high = zeros( NNk, 1 );

v_t_PF_c_low = zeros( NNk, 1 );
v_t_PF_c_med = zeros( NNk, 1 );
v_t_PF_c_high = zeros( NNk, 1 );

v_t_ADPTV_c_low = zeros( NNk, 1 );
v_t_ADPTV_c_med = zeros( NNk, 1 );
v_t_ADPTV_c_high = zeros( NNk, 1 );

v_t_PFmp_c_low = zeros( NNk, 1 );
v_t_PFmp_c_med = zeros( NNk, 1 );
v_t_PFmp_c_high = zeros( NNk, 1 );

v_t_DIM_low = zeros( NNk, 1 );
v_t_DIM_med = zeros( NNk, 1 );
v_t_DIM_high = zeros( NNk, 1 );




ca_v_c_CL_mp_ref = cell( NNk, 1  );
ca_v_c_CL_low = cell( NNk, 1  );
ca_v_c_CL_med = cell( NNk, 1  );
ca_v_c_CL_high = cell( NNk, 1  );

ca_v_c_SH_low = cell( NNk, 1  );
ca_v_c_SH_high = cell( NNk, 1  );

ca_v_c_PF_interp_low = cell( NNk, 1  );
ca_v_c_PF_interp_med = cell( NNk, 1  );
ca_v_c_PF_interp_high = cell( NNk, 1  );

ca_v_c_ADPTV_low = cell( NNk, 1  );
ca_v_c_ADPTV_med = cell( NNk, 1  );
ca_v_c_ADPTV_high = cell( NNk, 1  );

ca_v_c_PFmp_interp_low = cell( NNk, 1  );
ca_v_c_PFmp_interp_med = cell( NNk, 1  );
ca_v_c_PFmp_interp_high = cell( NNk, 1  );

ca_v_tc_dimM_low = cell( NNk, 1  );
ca_v_tc_dimM_med = cell( NNk, 1  );
ca_v_tc_dimM_high = cell( NNk, 1  );


for lp_Nk=1:NNk
       
    ca_v_c_CL_mp_ref{ lp_Nk } = ca_st_ref_results{ lp_Nk }.v_c_CL_mp_ref;

    v_t_CL_c_low( lp_Nk ) = ca_st_results{ lp_Nk }.t_CL_c_low;
    v_t_CL_c_med( lp_Nk ) = ca_st_results{ lp_Nk }.t_CL_c_med;
    v_t_CL_c_high( lp_Nk ) = ca_st_results{ lp_Nk }.t_CL_c_high;
    
    v_t_SH_c_low( lp_Nk ) = ca_st_results{ lp_Nk }.t_SH_c_low;
    v_t_SH_c_high( lp_Nk ) = ca_st_results{ lp_Nk }.t_SH_c_high;
    
    v_t_PF_c_low( lp_Nk ) = ca_st_results{ lp_Nk }.t_PF_c_low;
    v_t_PF_c_med( lp_Nk ) = ca_st_results{ lp_Nk }.t_PF_c_med;
    v_t_PF_c_high( lp_Nk ) = ca_st_results{ lp_Nk }.t_PF_c_high;
    
    v_t_ADPTV_c_low( lp_Nk ) = ca_st_results{ lp_Nk }.t_ADPTV_c_low;
    v_t_ADPTV_c_med( lp_Nk ) = ca_st_results{ lp_Nk }.t_ADPTV_c_med;
    v_t_ADPTV_c_high( lp_Nk ) = ca_st_results{ lp_Nk }.t_ADPTV_c_high;
    
    v_t_PFmp_c_low( lp_Nk ) = ca_st_results{ lp_Nk }.t_PFmp_c_low;
    v_t_PFmp_c_med( lp_Nk ) = ca_st_results{ lp_Nk }.t_PFmp_c_med;
    v_t_PFmp_c_high( lp_Nk ) = ca_st_results{ lp_Nk }.t_PFmp_c_high;
       
    v_t_DIM_low( lp_Nk ) = ca_st_results{ lp_Nk }.t_DIM_low;
    v_t_DIM_med( lp_Nk ) = ca_st_results{ lp_Nk }.t_DIM_med;
    v_t_DIM_high( lp_Nk ) = ca_st_results{ lp_Nk }.t_DIM_high;
    
    
    
    ca_v_c_CL_low{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_CL_low;
    ca_v_c_CL_med{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_CL_med;
    ca_v_c_CL_high{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_CL_high;
    
    ca_v_c_SH_low{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_SH_low;
    ca_v_c_SH_high{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_SH_high;
    
    ca_v_c_PF_interp_low{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_PF_interp_low;
    ca_v_c_PF_interp_med{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_PF_interp_med;
    ca_v_c_PF_interp_high{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_PF_interp_high;
    
    ca_v_c_ADPTV_low{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_ADPTV_low;
    ca_v_c_ADPTV_med{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_ADPTV_med;
    ca_v_c_ADPTV_high{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_ADPTV_high;

    ca_v_c_PFmp_interp_low{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_PFmp_interp_low;
    ca_v_c_PFmp_interp_med{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_PFmp_interp_med;
    ca_v_c_PFmp_interp_high{ lp_Nk } = ca_st_results{ lp_Nk }.v_c_PFmp_interp_high;
    
    ca_v_tc_dimM_low{ lp_Nk } = ca_st_results{ lp_Nk }.v_tc_dimM_low;
    ca_v_tc_dimM_med{ lp_Nk } = ca_st_results{ lp_Nk }.v_tc_dimM_med;
    ca_v_tc_dimM_high{ lp_Nk } = ca_st_results{ lp_Nk }.v_tc_dimM_high;
    
end


st_perf_tests_c = struct;
st_perf_tests_c.v_Nk = v_Nk;
st_perf_tests_c.ca_v_k_q = ca_v_k_q;
st_perf_tests_c.k_min = k_min;
st_perf_tests_c.k_max = k_max;
st_perf_tests_c.st_r_shear = st_r_shear;
st_perf_tests_c.st_p = st_p;

st_perf_tests_c.v_t_CL_c_low = v_t_CL_c_low;
st_perf_tests_c.v_t_CL_c_med = v_t_CL_c_med;
st_perf_tests_c.v_t_CL_c_high = v_t_CL_c_high;

st_perf_tests_c.v_t_SH_c_low = v_t_SH_c_low;
st_perf_tests_c.v_t_SH_c_high = v_t_SH_c_high;

st_perf_tests_c.v_t_PF_c_low = v_t_PF_c_low;
st_perf_tests_c.v_t_PF_c_med = v_t_PF_c_med;
st_perf_tests_c.v_t_PF_c_high = v_t_PF_c_high;

st_perf_tests_c.v_t_ADPTV_c_low = v_t_ADPTV_c_low;
st_perf_tests_c.v_t_ADPTV_c_med = v_t_ADPTV_c_med;
st_perf_tests_c.v_t_ADPTV_c_high = v_t_ADPTV_c_high;

st_perf_tests_c.v_t_PFmp_c_low = v_t_PFmp_c_low;
st_perf_tests_c.v_t_PFmp_c_med = v_t_PFmp_c_med;
st_perf_tests_c.v_t_PFmp_c_high = v_t_PFmp_c_high;

st_perf_tests_c.v_t_DIM_low = v_t_DIM_low;
st_perf_tests_c.v_t_DIM_med = v_t_DIM_med;
st_perf_tests_c.v_t_DIM_high = v_t_DIM_high;



st_perf_tests_c.ca_v_c_CL_mp_ref = ca_v_c_CL_mp_ref;
st_perf_tests_c.ca_v_c_CL_low = ca_v_c_CL_low;
st_perf_tests_c.ca_v_c_CL_med = ca_v_c_CL_med;
st_perf_tests_c.ca_v_c_CL_high = ca_v_c_CL_high;

st_perf_tests_c.ca_v_c_SH_low = ca_v_c_SH_low;
st_perf_tests_c.ca_v_c_SH_high = ca_v_c_SH_high;

st_perf_tests_c.ca_v_c_PF_interp_low = ca_v_c_PF_interp_low;
st_perf_tests_c.ca_v_c_PF_interp_med = ca_v_c_PF_interp_med;
st_perf_tests_c.ca_v_c_PF_interp_high = ca_v_c_PF_interp_high;

st_perf_tests_c.ca_v_c_ADPTV_low = ca_v_c_ADPTV_low;
st_perf_tests_c.ca_v_c_ADPTV_med = ca_v_c_ADPTV_med;
st_perf_tests_c.ca_v_c_ADPTV_high = ca_v_c_ADPTV_high;

st_perf_tests_c.ca_v_c_PFmp_interp_low = ca_v_c_PFmp_interp_low;
st_perf_tests_c.ca_v_c_PFmp_interp_med = ca_v_c_PFmp_interp_med;
st_perf_tests_c.ca_v_c_PFmp_interp_high = ca_v_c_PFmp_interp_high;

st_perf_tests_c.ca_v_tc_dimM_low = ca_v_tc_dimM_low;
st_perf_tests_c.ca_v_tc_dimM_med = ca_v_tc_dimM_med;
st_perf_tests_c.ca_v_tc_dimM_high = ca_v_tc_dimM_high;


save( s_filename, 'st_perf_tests_c' );








end
function fn_ww__sim_re__acc_tests_c__red( Nz_start, Nz_end, b_calc_ref )
%fn_ww__sim_re__acc_tests_c__red: Simulation CL+PF eig accuracy, backwards error, cond number tests
%
%   fn_ww__sim_re__acc_tests_c__red( Nz_start, Nz_end, b_calc_ref )
%
% Calculate accuracy for eigenvalues from collocation and path-following
% calcs of Rayleigh equation.
%
% TAGS: SISCPFLIB
%
% See also
%   fn_ww__sim_re__acc_tests_c__red__parfor(),
%   fn_ww__analyse_re__acc_tests_c__red()


[ st_p ] = fn_ww__setup__param_std__re_cl(  );
phy_h = 20; Fr2 = 0.05;
[ st_fn_shear, st_p ] = fn_ww__setup__shear_fn__nondim_cospwr( st_p, sqrt( Fr2 * st_p.phy_g * phy_h ), phy_h, 2, 4 * pi, 1, 0.5 );
[ st_r_shear ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear, st_p );

% Nz_ref = 384;
% mp_digits_ref = 56;
Nz_ref = 64;
mp_digits_ref = 34;
k_min = 0.025;
k_max = 250;
Nk = 12;

NNz = Nz_end - Nz_start + 1;
v_Nz = Nz_start:Nz_end;




if ( b_calc_ref )
    
    % Calculate reference data
    
    % Create v_k vector, log distributed
    v_k = fn_ww__util__create_k_vec( k_min, k_max, Nk, 1, 0 );    
    
    % Setup matrices
    [ st_Dn_mp_ref ] = fn_ww__setup__diffmtrx_mp__WR_chebdif( Nz_ref, mp_digits_ref );
    [ st_Dn_mp_ref ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn_mp_ref, st_p );    
    %[ st_v_shear_mp_ref ] = fn_ww__setup__shear_fn_to_vec( st_Dn_mp_ref, st_fn_shear, st_p );
    
    % Do calc    
    [ st_p_ref ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_mp', true, 'bp_disp_update', true, 'ip_mp_digits', mp_digits_ref ) );
    fprintf( 'Starting MP ref calc... \n' );
    %[ v_c_CL_ref ] = fn_ww__calc_re__cl__red_c( st_Dn_mp_ref, v_k, st_v_shear_mp_ref, st_p_ref );
    [ v_c_CL_ref ] = fn_ww__calc_re__cl__red_c( st_Dn_mp_ref, v_k, st_r_shear, st_p_ref );
    fprintf( '... end.\n' );
    
    % Save
    st_ref = struct;
    st_ref.v_k = v_k;
    st_ref.v_c_CL_ref = v_c_CL_ref;
    st_ref.st_Dn_mp_ref = st_Dn_mp_ref;
    st_ref.st_v_shear_mp_ref = st_v_shear_mp_ref;
    st_ref.st_p_ref = st_p_ref;
    st_ref.Nk = Nk;
    st_ref.k_min = k_min;
    st_ref.k_max = k_max;
    st_ref.Nz_ref = Nz_ref;
    
    save( 'st_acc_tests__ref', 'st_ref' );   
    
else

    % TODO should actually load this
    
    % Create v_k vector, log distributed
    v_k = fn_ww__util__create_k_vec( k_min, k_max, Nk, 1, 0 );
   
end

v_k_phy = v_k / st_p.phy_h;

[ st_adptv ] = fn_ww__calc_re__adptv_interval_refine( v_k, st_p );  



ca_v_c_CL = cell( NNz, 1 );
ca_v_c_PF = cell( NNz, 1 );
ca_v_c_PFmp = cell( NNz, 1 );
ca_v_c_ADPTV_CL = cell( NNz, 1 );
ca_v_c_ADPTV_PF = cell( NNz, 1 );
ca_v_c_ADPTV_PFmp = cell( NNz, 1 );
ca_v_c_DIM = cell( NNz, 1 );
ca_st_err_CL = cell( NNz, 1 );
ca_st_err_PF = cell( NNz, 1 );

parpool(4);
parfor lp_Nz=1:NNz

    fprintf( '--[ Starting %d ]--\n', lp_Nz );
    
    % Do the calculation
    [ ca_v_c_CL{lp_Nz,1}, ...
        ca_v_c_PF{lp_Nz,1}, ...
        ca_v_c_PFmp{lp_Nz,1}, ...
        ca_v_c_ADPTV_CL{lp_Nz,1}, ...
        ca_v_c_ADPTV_PF{lp_Nz,1}, ...
        ca_v_c_ADPTV_PFmp{lp_Nz,1}, ...        
        ca_v_c_DIM{lp_Nz,1}, ...
        ca_st_err_CL{lp_Nz,1}, ...
        ca_st_err_PF{lp_Nz,1} ] = fn_ww__sim_re__acc_tests_c__red__parfor( v_Nz(lp_Nz), st_r_shear, st_adptv, k_min, k_max, v_k, v_k_phy, st_p );
    
    fprintf( '--[ End %d ]--\n', lp_Nz );
    
end
delete(gcp);

st_data = struct;
st_data.ca_v_c_CL = ca_v_c_CL;
st_data.ca_v_c_PF = ca_v_c_PF;
st_data.ca_v_c_PFmp = ca_v_c_PFmp;
st_data.ca_v_c_DIM = ca_v_c_DIM;
st_data.ca_st_err_CL = ca_st_err_CL;
st_data.ca_st_err_PF = ca_st_err_PF;
st_data.v_k = v_k;
st_data.v_k_phy = v_k_phy;
st_data.st_p = st_p;
st_data.NNz = NNz;
st_data.v_Nz = v_Nz;
st_data.st_fn_shear = st_fn_shear;
st_data.k_min = k_min;
st_data.k_max = k_max;

save( 'st_acc_tests__data', 'st_data' );




end
function [ st_results ] = fn_ww__sim_re__perf_tests_c__red__parfor( st_data, v_k_q, Nrep, b_do_ref, b_do_CL, b_do_DIMmed, b_do_DIMhigh, b_no_save, st_p )
%fn_ww__sim_re__perf_tests_c__red__parfor: Simulation (parfor fn) for CL,DIM,PF perf tests
%
%   [ st_results ] = fn_ww__sim_re__perf_tests_c__red__parfor( st_data, Nk, Nrep, b_do_ref, b_do_CL, b_do_DIMmed, b_do_DIMhigh, b_no_save, st_p )
%
% Only call from fn_ww__sim_re__perf_tests_c__red(). Performs the actual
% test in a separate function so that it can be executed within a parfor
% loop.
% 
% TAGS: SISCPFLIB
%
% INPUT
% 
% st_data : struct, containing data to setup shear profiles, etc.
% 
% Nk : Number of k elements to use
% 
% Nrep : Repetitions for each algorithm's test
% 
% b_do_CL : Flag, do the collocation tests
% 
% b_no_save : Flag, whether to not save the actual c vectors
% 
% OUTPUT
% 
% st_results : struct, containing all the necessary
%
%
% See also
%   FN_WW__SIM_RE__PERF_TESTS_C__RED(), 
%   FN_WW__SIM_RE__PERF_TESTS_C__RED__SET()
% 
% 






    % Setup
    Nk = numel( v_k_q );
    v_k_q_phy = v_k_q / st_p.phy_h;
    PF_epsilon = 1e-5;
    st_results = struct;
    [ st_adptv ] = fn_ww__calc_re__adptv_interval_refine( v_k_q, st_p );    

    % Make things simplier... assume CL doesn't need reps over threshold
    NCLrep = Nrep;
    if ( Nk > 500 )
        NCLrep = min( [ 3 Nrep ] );
    end    
    
    
    % Debug test (MATLAB doing weird stuff)
    t_PF_c_low = 0;
    t_PF_c_med = 0;
    t_PF_c_high = 0;
    st_results.t_PF_c_low = 0;
    st_results.t_PF_c_med= 0;
    st_results.t_PF_c_high = 0;
    
    t_ADPTV_c_low = 0;
    t_ADPTV_c_med = 0;
    t_ADPTV_c_high = 0;
    st_results.t_ADPTV_c_low = 0;
    st_results.t_ADPTV_c_med = 0;
    st_results.t_ADPTV_c_high = 0;
    
    t_PFmp_c_low = 0;
    t_PFmp_c_med = 0;
    t_PFmp_c_high = 0;
    st_results.t_PFmp_c_low = 0;
    st_results.t_PFmp_c_med = 0;
    st_results.t_PFmp_c_high = 0;     
    
    t_CL_c_low = 0;
    t_CL_c_med = 0;
    t_CL_c_high = 0;
    st_results.t_CL_c_low = 0;
    st_results.t_CL_c_med = 0;
    st_results.t_CL_c_higjh = 0;
    
    t_SH_c_low = 0;
    t_SH_c_high = 0;    
    st_results.t_SH_c_low = 0;
    st_results.t_SH_c_high = 0;
    
    t_DIM_low = 0;
    t_DIM_med = 0;
    t_DIM_high = 0;
    st_results.t_DIM_low = 0;
    st_results.t_DIM_med = 0;
    st_results.t_DIM_high = 0;
    
    
    
    
    
    %--[ PF tests ]--------------------------------------------------------    
    % Low
    t_vlarge = 1e50;
    for lp_r=1:Nrep
        tic
        [ st_dp_PF_low ] = fn_ww__calc_re__clpf_dp__red_c( st_data.st_Dn_std_low, st_data.k_min-PF_epsilon, st_data.k_max+PF_epsilon, st_data.st_r_shear, 2.0, 1e-5, st_p );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end 
    t_PF_c_low = t_vlarge;
    t_vlarge = 1e50;    
    for lp_r=1:Nrep
        tic
        [ v_c_PF_interp_low, ~ ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_PF_low, v_k_q, 0 );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end    
    st_results.t_PF_c_low = t_PF_c_low + t_vlarge;
 
    % If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear st_dp_PF_low;
        clear v_c_PF_interp_low;
        st_dp_PF_low = -1;
        v_c_PF_interp_low = -1;
    end
    
    % N=48, tol=11
    t_vlarge = 1e50;
    for lp_r=1:Nrep
        tic
        [ st_dp_PF_med ] = fn_ww__calc_re__clpf_dp__red_c( st_data.st_Dn_std_med, st_data.k_min-PF_epsilon, st_data.k_max+PF_epsilon, st_data.st_r_shear, 2.0, 1e-8, st_p );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end 
    t_PF_c_med = t_vlarge;
    t_vlarge = 1e50;    
    for lp_r=1:Nrep
        tic
        [ v_c_PF_interp_med, ~ ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_PF_med, v_k_q, 0 );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end    
    st_results.t_PF_c_med = t_PF_c_med + t_vlarge;
    
    % If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear st_dp_PF_med;
        clear v_c_PF_interp_med;
        st_dp_PF_med = -1;
        v_c_PF_interp_med = -1;
    end    
    
    % N=64, tol=13
    t_vlarge = 1e50;
    for lp_r=1:Nrep
        tic
        [ st_dp_PF_high ] = fn_ww__calc_re__clpf_dp__red_c( st_data.st_Dn_std_high, st_data.k_min-PF_epsilon, st_data.k_max+PF_epsilon, st_data.st_r_shear, 2.0, 1e-11, st_p );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end 
    t_PF_c_high = t_vlarge;
    t_vlarge = 1e50;    
    for lp_r=1:Nrep
        tic
        [ v_c_PF_interp_high, ~ ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_PF_high, v_k_q, 0 );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end    
    st_results.t_PF_c_high = t_PF_c_high + t_vlarge;
    
    % If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear st_dp_PF_high;
        clear v_c_PF_interp_high;
        st_dp_PF_high = -1;
        v_c_PF_interp_high = -1;
    end    
    

    
    
    
    %--[ Adptv PF ]--------------------------------------------------------
    % N=16
    t_vlarge = 1e50;
    st_proc = struct;
    st_proc.ip_method = 2;
    st_proc.pf_tol = 1e-5;
    for lp_r=1:Nrep
        tic
        [ v_c_ADPTV_low ] = fn_ww__calc_re__adptv__red_c( st_data.st_Dn_adptv_low, st_data.st_r_shear, st_adptv, st_proc, st_p );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end
    st_results.t_ADPTV_c_low = t_vlarge;
                 
	% If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear v_c_ADPTV_low;
        v_c_ADPTV_low = -1;
    end
    
    % N=32
    t_vlarge = 1e50;
    st_proc = struct;
    st_proc.ip_method = 2;
    st_proc.pf_tol = 1e-8;
    for lp_r=1:Nrep
        tic
        [ v_c_ADPTV_med ] = fn_ww__calc_re__adptv__red_c( st_data.st_Dn_adptv_med, st_data.st_r_shear, st_adptv, st_proc, st_p );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end
    st_results.t_ADPTV_c_med = t_vlarge;
                 
	% If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear v_c_ADPTV_med;
        v_c_ADPTV_med = -1;
    end
    
    % N=48
    t_vlarge = 1e50;
    st_proc = struct;
    st_proc.ip_method = 2;
    st_proc.pf_tol = 1e-11;
    for lp_r=1:Nrep
        tic
        [ v_c_ADPTV_high ] = fn_ww__calc_re__adptv__red_c( st_data.st_Dn_adptv_high, st_data.st_r_shear, st_adptv, st_proc, st_p );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end
    st_results.t_ADPTV_c_high = t_vlarge;
                 
	% If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear v_c_ADPTV_high;
        v_c_ADPTV_high = -1;
    end
    
    
    
    
    
    
    %--[ PFmp tests ]--------------------------------------------------                
    % 
    t_vlarge = 1e50;
    for lp_r=1:Nrep
        tic
        [ st_p_PFmp_low ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_pf_use_eig_mp', true, 'st_Dn_mp', st_data.st_Dn_mp_low ) );
        [ st_dp_PFmp_low ] = fn_ww__calc_re__clpf_dp__red_c( st_data.st_Dn_std_low, st_data.k_min-PF_epsilon, st_data.k_max+PF_epsilon, st_data.st_r_shear, 2.0, 1e-5, st_p_PFmp_low );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end 
    t_PFmp_c_low = t_vlarge;
    t_vlarge = 1e50;    
    for lp_r=1:Nrep
        tic
        [ v_c_PFmp_interp_low, ~ ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_PFmp_low, v_k_q, 0 );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end    
    st_results.t_PFmp_c_low = t_PFmp_c_low + t_vlarge;
 
    % If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear st_dp_PFmp_low;
        clear v_c_PFmp_interp_low;
        st_dp_PFmp_low = -1;
        v_c_PFmp_interp_low = -1;
    end
    
    % N=48, tol=11   
    t_vlarge = 1e50;
    for lp_r=1:Nrep
        tic
        [ st_p_PFmp_med ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_pf_use_eig_mp', true, 'st_Dn_mp', st_data.st_Dn_mp_med ) );
        [ st_dp_PFmp_med ] = fn_ww__calc_re__clpf_dp__red_c( st_data.st_Dn_std_med, st_data.k_min-PF_epsilon, st_data.k_max+PF_epsilon, st_data.st_r_shear, 2.0, 1e-8, st_p_PFmp_med );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end 
    t_PFmp_c_med = t_vlarge;
    t_vlarge = 1e50;    
    for lp_r=1:Nrep
        tic
        [ v_c_PFmp_interp_med, ~ ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_PFmp_med, v_k_q, 0 );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end    
    st_results.t_PFmp_c_med = t_PFmp_c_med + t_vlarge;
    
    % If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear st_dp_PFmp_med;
        clear v_c_PFmp_interp_med;
        st_dp_PFmp_med = -1;
        v_c_PFmp_interp_med = -1;
    end    
    
    % N=64, tol=13
    t_vlarge = 1e50;
    for lp_r=1:Nrep
        tic
        [ st_p_PF64_high ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_pf_use_eig_mp', true, 'st_Dn_mp', st_data.st_Dn_mp_high ) );
        [ st_dp_PFmp_high ] = fn_ww__calc_re__clpf_dp__red_c( st_data.st_Dn_std_high, st_data.k_min-PF_epsilon, st_data.k_max+PF_epsilon, st_data.st_r_shear, 2.0, 1e-11, st_p_PF64_high );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end 
    t_PFmp_c_high = t_vlarge;
    t_vlarge = 1e50;    
    for lp_r=1:Nrep
        tic
        [ v_c_PFmp_interp_high, ~ ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_PFmp_high, v_k_q, 0 );
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end    
    st_results.t_PFmp_high = t_PFmp_c_high + t_vlarge;
    
    % If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear st_dp_PFmp_high;
        clear v_c_PFmp_interp_high;
        st_dp_PFmp_high = -1;
        v_c_PFmp_interp_high = -1;
    end    
    
    
    
    
    
    
    
    
  
    
    
    
    
    %--[ CL ]--------------------------------------------------------------
    if ( b_do_CL )

        [ st_p_CL ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_disp_update', false, 'bp_mp', false ) );
           
        % N=32
        t_vlarge = 1e50;
        st_proc = struct;
        st_proc.ip_method = 1;
        for lp_r=1:NCLrep
            tic
            [ v_c_CL_low ] = fn_ww__calc_re__adptv__red_c( st_data.st_Dn_adptv_low, st_data.st_r_shear, st_adptv, st_proc, st_p );
            %[ v_c_CL_low ] = fn_ww__calc_re__cl__red_c( st_data.st_Dn_std_low, v_k_q, st_data.st_r_shear, st_p_CL );
            t_new = toc;
            if ( t_new < t_vlarge ), t_vlarge = t_new; end
        end
        st_results.t_CL_c_low = t_vlarge;

        % N=48
        t_vlarge = 1e50;
        st_proc = struct;
        st_proc.ip_method = 1;
        for lp_r=1:NCLrep
            tic
            [ v_c_CL_med ] = fn_ww__calc_re__adptv__red_c( st_data.st_Dn_adptv_med, st_data.st_r_shear, st_adptv, st_proc, st_p );
            %[ v_c_CL_med ] = fn_ww__calc_re__cl__red_c( st_data.st_Dn_std_med, v_k_q, st_data.st_r_shear, st_p_CL );
            t_new = toc;
            if ( t_new < t_vlarge ), t_vlarge = t_new; end
        end
        st_results.t_CL_c_med = t_vlarge;

        % N=64
        t_vlarge = 1e50;
        st_proc = struct;
        st_proc.ip_method = 1;
        for lp_r=1:NCLrep
            tic
            [ v_c_CL_high ] = fn_ww__calc_re__adptv__red_c( st_data.st_Dn_adptv_high, st_data.st_r_shear, st_adptv, st_proc, st_p );
            %[ v_c_CL_high ] = fn_ww__calc_re__cl__red_c( st_data.st_Dn_std_high, v_k_q, st_data.st_r_shear, st_p_CL );
            t_new = toc;
            if ( t_new < t_vlarge ), t_vlarge = t_new; end
        end    
        st_results.t_CL_c_high = t_vlarge;    
                
    else
        
        % Since we haven't tested, set to a value we can filter out later.
        v_c_CL_low = -1;
        v_c_CL_med = -1;
        v_c_CL_high = -1;        
        st_results.t_CL_c_low = -1;
        st_results.t_CL_c_med = -1;
        st_results.t_CL_c_high = -1;
        
    end
    

    % If we don't need to save the data (updated this 20180815, hoping I'm
    % doing what I originally intended...)
    if ( b_no_save ~= 0 )
        clear v_c_CL_low;
        clear v_c_CL_med;
        clear v_c_CL_high;
        
        v_c_CL_low = -1;
        v_c_CL_med = -1;
        v_c_CL_high = -1;
    end        
    
    


    
    
    %--[ DIM tests ]-------------------------------------------------------
    % LOW
    t_vlarge = 1e50;
    for lp_r=1:Nrep   
        tic
        [ ~, ~, v_tc_dimM_low, ~, ~ ] = fn_ww__ext__dim__c( v_k_q_phy, 0, st_data.st_fn_shear_DIM, 1e-10, 28, 10, st_p.phy_h, st_p.phy_g, false );
        v_tc_dimM_low = v_tc_dimM_low / st_p.phy_U0;
        t_new = toc;
        if ( t_new < t_vlarge ), t_vlarge = t_new; end
    end    
    st_results.t_DIM_low = t_vlarge;

    % If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear v_tc_dimM_low;
        v_tc_dimM_low = -1;
    end    
    
    
    if ( b_do_DIMmed )
        % MED
        t_vlarge = 1e50;
        for lp_r=1:Nrep    
            tic
            [ ~, ~, v_tc_dimM_med, ~, ~ ] = fn_ww__ext__dim__c( v_k_q_phy, 0, st_data.st_fn_shear_DIM, 1e-12, 800, 10, st_p.phy_h, st_p.phy_g, false );
            v_tc_dimM_med = v_tc_dimM_med / st_p.phy_U0;
            t_new = toc;
            if ( t_new < t_vlarge ), t_vlarge = t_new; end
        end    
        st_results.t_DIM_med = t_vlarge;
    else        
        st_results.t_DIM_med = -1;
        v_tc_dimM_med = -1;
    end
    
    % If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear v_tc_dimM_med;
        v_tc_dimM_med = -1;
    end
    
    
    if ( b_do_DIMhigh )
        % HIGH
        t_vlarge = 1e50;
        for lp_r=1:Nrep   
            tic
            [ ~, ~, v_tc_dimM_high, ~, ~ ] = fn_ww__ext__dim__c( v_k_q_phy, 0, st_data.st_fn_shear_DIM, 1e-13, 23000, 20, st_p.phy_h, st_p.phy_g, false );
            v_tc_dimM_high = v_tc_dimM_high / st_p.phy_U0;
            t_new = toc;
            if ( t_new < t_vlarge ), t_vlarge = t_new; end
        end    
        st_results.t_DIM_high = t_vlarge;
    else        
        st_results.t_DIM_high = -1;
        v_tc_dimM_high = -1;        
    end
        
    % If we don't need to save the data, the memory
    if ( b_no_save ~= 0 )
        clear v_tc_dimM_high;
        v_tc_dimM_high = -1;
    end    
    

    
    
    

    %--[ DIM tests ]-------------------------------------------------------
    % Just use flag for DIMhigh
    if ( b_do_DIMhigh )   
    
        % Shooting med
        t_vlarge = 1e50;
        st_tol = struct;
        st_tol.tol_fzero = 1e-6;
        st_tol.tol_ode45_rel = 5e-6;
        st_tol.tol_ode45_abs = 1e-6;
        st_tol.tol_depth = 1e-7;
        for lp_r=1:NCLrep
            tic
            [ v_c_SH_low ] = fn_ww__calc_shoot__red_c( v_k_q, st_data.st_r_shear, st_tol, st_p );
            t_new = toc;
            if ( t_new < t_vlarge ), t_vlarge = t_new; end
        end    
        st_results.t_SH_c_low = t_vlarge;

        % Shooting high
        t_vlarge = 1e50;
        st_tol = struct;
        st_tol.tol_fzero = 1e-10;
        st_tol.tol_ode45_rel = 1e-8;
        st_tol.tol_ode45_abs = 1e-9;
        st_tol.tol_depth = 1e-10;
        for lp_r=1:NCLrep
            tic
            [ v_c_SH_high ] = fn_ww__calc_shoot__red_c( v_k_q, st_data.st_r_shear, st_tol, st_p );
            t_new = toc;
            if ( t_new < t_vlarge ), t_vlarge = t_new; end
        end    
        st_results.t_SH_c_high = t_vlarge;
    
    else
        
        % Since we haven't tested, set to a value we can filter out later.
        v_c_SH_low = -1;
        v_c_SH_high = -1;
        st_results.t_SH_c_low = -1;
        st_results.t_SH_c_high = -1;
    
    end    
    
    if ( b_no_save ~= 0 )
        clear v_c_SH_low;
        clear v_c_SH_high;
        
        v_c_SH_low = -1;
        v_c_SH_high = -1;
    end         
    
    
    
    
    
    
    %--[ Store results ]---------------------------------------------------
    st_results.v_k_q_phy = v_k_q_phy;    
   
    st_results.v_c_CL_low = v_c_CL_low;
    st_results.v_c_CL_med = v_c_CL_med;
    st_results.v_c_CL_high = v_c_CL_high;
    
    st_results.v_c_SH_low = v_c_SH_low;
    st_results.v_c_SH_high = v_c_SH_high;
    
    st_results.v_c_PF_interp_low = v_c_PF_interp_low;
    st_results.v_c_PF_interp_med = v_c_PF_interp_med; 
    st_results.v_c_PF_interp_high = v_c_PF_interp_high;

    st_results.v_c_ADPTV_low = v_c_ADPTV_low;
    st_results.v_c_ADPTV_med = v_c_ADPTV_med;
    st_results.v_c_ADPTV_high = v_c_ADPTV_high;
    
    st_results.v_c_PFmp_interp_low = v_c_PFmp_interp_low;
    st_results.v_c_PFmp_interp_med = v_c_PFmp_interp_med; 
    st_results.v_c_PFmp_interp_high = v_c_PFmp_interp_high;    
    
    st_results.v_tc_dimM_low = v_tc_dimM_low;
    st_results.v_tc_dimM_med = v_tc_dimM_med;
    st_results.v_tc_dimM_high = v_tc_dimM_high;  
    

end
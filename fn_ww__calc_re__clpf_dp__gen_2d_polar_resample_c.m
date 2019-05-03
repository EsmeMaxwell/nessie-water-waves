function [ st_gen_resample ] = fn_ww__calc_re__clpf_dp__gen_2d_polar_resample_c( st_Dn, k_min, k_max, Nkresample, Ntheta, st_v_shear_x, st_v_shear_y, k0, tol, st_p )
%fn_ww__calc_re__clpf_dp__gen_2d_polar_resample_c: INCOMPLETE Calc PF-c disperion relation for general scattered data
% 
%   [ st_gen_resample ] = fn_ww__calc_re__clpf_dp__gen_2d_polar_resample_c( st_Dn, k_min, k_max, Nkresample, Ntheta, st_v_shear_x, st_v_shear_y, k0, tol, st_p )
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%



a_D2m = st_Dn.a_D2m;
a_Dm = st_Dn.a_Dm;

N = numel( st_Dn.v_z0 );
v_theta = fn_ww__util__create_theta_vec( Ntheta );
v_k = fn_ww__util__create_k_vec( k_min, k_max, Nkresample, 3, 1e-4 );

a_I = eye( size( a_Dm ) );
a_0I = 0 * a_I;
a_R = a_I;
a_R(1,1) = 0;
a_R(end,end) = 0;
a_R_short = a_R(1:end-1,1:end-1);

% First step, create angular DP result
[ st_dp_angular ] = fn_ww__calc_re__clpf_dp__gen_sl_polar_ang_c( st_Dn, k0, st_v_shear_x, st_v_shear_y, tol, st_p );

% Generate the v_y0 starting vector in angular circle, used for initial y0
% in radial slice calculations
a_y0_ang_PF_interp = zeros( N+1, Ntheta );
[ a_y0_ang_PF_interp ] = fn_ww__calc_re__clpf_interp__sl_ang_c( st_dp_angular, v_theta, 1 );

% Init
a_c_slice_radial_resampled = zeros( Ntheta, Nkresample );
a_dk_c_slice_radial_resampled = zeros( Ntheta, Nkresample );
a_dtheta_c_slice_radial_resampled = zeros( Ntheta, Nkresample );

ca_st_dp_radial = cell( 1, Ntheta );
for lp_theta=1:Ntheta
    
    theta = v_theta( lp_theta );
    
%    fprintf( '%d of %d, theta=%f\n', lp_theta, Ntheta, theta );
    
    % Calculate DP control points along radial slice
    [ st_p_for_rad_c ] = fn_ww__setup__merge_parameters( st_p, struct( 'fp_k0', k0, 'vp_pf_v_y0', a_y0_ang_PF_interp( :, lp_theta ) ) );   
    [ st_dp_radial ] = fn_ww__calc_re__clpf_dp__gen_sl_polar_rad_c( st_Dn, k_min, k_max, theta, st_v_shear_x, st_v_shear_y, tol, st_p_for_rad_c );
    
    % Save it, we use it in debuging (my omit this if wanting optimisation)
    ca_st_dp_radial{lp_theta} = st_dp_radial;
    
    % Calculate resampled k points using the interpolant (we need the
    % vector form to be able to calculate the angular derivative further
    % below)    
    [ a_y0_slice_radial_resampled, a_dy0_slice_radial_resampled ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp_radial, v_k, 1 );
    a_c_slice_radial_resampled(lp_theta,:) = a_y0_slice_radial_resampled( end, : );
    a_dk_c_slice_radial_resampled(lp_theta,:) = a_dy0_slice_radial_resampled( end, : );  
    
    % Calculate angular derivative too
    for lp_k=1:Nkresample
        [ ~, a_dtheta_c_slice_radial_resampled( lp_theta, lp_k ) ] = fn_local__calc_dtheta( v_k(lp_k), theta, a_y0_slice_radial_resampled( 1:end-1, lp_k ), a_y0_slice_radial_resampled( end, lp_k ) );
    end
            
end



st_gen_resample = struct;
st_gen_resample.st_dp_angular = st_dp_angular;
st_gen_resample.a_y0_ang_PF_interp = a_y0_ang_PF_interp;
st_gen_resample.ca_st_dp_radial = ca_st_dp_radial;
st_gen_resample.a_c_slice_radial_resampled = a_c_slice_radial_resampled;
st_gen_resample.a_dk_c_slice_radial_resampled = a_dk_c_slice_radial_resampled;
st_gen_resample.a_dtheta_c_slice_radial_resampled = a_dtheta_c_slice_radial_resampled;
st_gen_resample.v_theta = v_theta;
st_gen_resample.v_k = v_k;
st_gen_resample.k_min = k_min;
st_gen_resample.k_max = k_max;
st_gen_resample.Nkresample = Nkresample;
st_gen_resample.Ntheta = Ntheta;
st_gen_resample.st_v_shear_x = st_v_shear_x;
st_gen_resample.st_v_shear_y = st_v_shear_y;
st_gen_resample.k0 = k0;
% st_gen_resample.a = a;
% st_gen_resample.b = b;
st_gen_resample.N = N;
st_gen_resanple.st_p = st_p;

return



    % Reused from DP angular code...
    function [ v_dtheta_w, dtheta_c ] = fn_local__calc_dtheta( k, theta, v_w, c )

        % U matrices
        v_U = cos(theta) * st_v_shear_x.v_U + sin(theta) * st_v_shear_y.v_U;
        v_dU = cos(theta) * st_v_shear_x.v_dU + sin(theta) * st_v_shear_y.v_dU;
        v_ddU = cos(theta) * st_v_shear_x.v_ddU + sin(theta) * st_v_shear_y.v_ddU;

        v_dtheta_U = -sin(theta) * st_v_shear_x.v_U + cos(theta) * st_v_shear_y.v_U;
        v_dtheta_dU = -sin(theta) * st_v_shear_x.v_dU + cos(theta) * st_v_shear_y.v_dU;
        v_dtheta_ddU = -sin(theta) * st_v_shear_x.v_ddU + cos(theta) * st_v_shear_y.v_ddU;

        a_U = diag( v_U );
        a_dU = diag( v_dU );
        a_ddU = diag( v_ddU );            

        a_dtheta_U = diag( v_dtheta_U );
        a_dtheta_dU = diag( v_dtheta_dU );
        a_dtheta_ddU = diag( v_dtheta_ddU );

        % LHS
        a_LHS = a_U * a_D2m - c * a_D2m - a_ddU - k^2 * a_U + c * k^2 * a_I;
        a_LHS_FS = a_U.^2 * a_Dm - 2 * c * a_U * a_Dm + c^2 * a_Dm - a_U * a_dU + c * a_dU - ( 1 / st_p.Fr2 ) * a_I;
        a_LHS(1,:) =  a_LHS_FS(1,:);

        % LHS extra
        a_LHSE = -a_D2m + k^2 * a_I;
        a_LHSE_FS = -2 * a_U * a_Dm + 2 * c * a_Dm + a_dU;
        a_LHSE(1,:) = a_LHSE_FS(1,:);

        % RHS
        a_RHS = -a_dtheta_U * a_D2m + a_dtheta_ddU + k^2 * a_dtheta_U;
        a_RHS_FS = -2 * a_U * a_dtheta_U * a_Dm + 2 * c * a_dtheta_U * a_Dm + a_dtheta_U * a_dU + a_U * a_dtheta_dU - c * a_dtheta_dU;
        a_RHS(1,:) = a_RHS_FS(1,:);

        % Shorten
        a_LHS = a_LHS(1:end-1,1:end-1);
        a_LHSE = a_LHSE(1:end-1,1:end-1);    
        a_RHS = a_RHS(1:end-1,1:end-1);
        %v_w = v_w(1:end-1);

        % Full M matrix
        a_M = [ a_LHS a_LHSE*v_w(1:end-1); -(a_R_short*v_w(1:end-1))' 0 ];
        a_L = [ a_RHS * v_w(1:end-1); 0 ];

        % Solve
        warning('');  %clear warning status
        v_dy = a_M \ a_L;        
        v_dtheta_w = v_dy(1:end-1);
        dtheta_c = v_dy(end);
    %     v_dy = [ v_dy; v_dy(end) ];  %fudge
    %     v_dy(end-1) = 0;
        [warnMsg, warnId ] = lastwarn;
        if ( ~isempty(warnMsg) )
            error( 'Problem ill-conditioned!');
        end

    end








end
function [ v_c_PF_interp ] = fn_ww__calc_re__clpf_interp__2d_resample_c( st_gen_resample, v_theta_q, v_k_q )
%fn_ww__calc_re__clpf_interp__2d_resample_c: INCOMPLETE Calc interpolant for PF-G scattered data
%
%   [ v_c_PF_interp ] = fn_ww__calc_re__clpf_interp__2d_resample_c( st_gen_resample, v_theta_q, v_k_q )
%
% v_theta_q, v_k_q : list of POINTS! (theta,k) style
%

if ( size( v_k_q, 1 ) ~= 1 )
    warning( 'Expected row vector for v_k_q' );
    v_k_q = v_k_q';
end

if ( size( v_theta_q, 1 ) ~= 1 )
    warning( 'Expected row vector for v_theta_q' );
    v_theta_q = v_theta_q';
end


v_theta = st_gen_resample.v_theta;
v_k = st_gen_resample.v_k;
k_min = st_gen_resample.k_min;
k_max = st_gen_resample.k_max;
Ntheta = st_gen_resample.Ntheta;
Nkresample = st_gen_resample.Nkresample;


% DO NOT need to worry about this any more, as we wrap back into the valid
% domain! :-)
% % Check theta vector is in correct format and append the -pi point onto pi
% % as the closing (with corresponding radial slice)
% assert( v_theta(1) == -pi, 'Start of v_theta range must be -pi.' );  % And YES we really do mean equality here.
% assert( v_theta(end) > -pi && v_theta(end) < pi, 'End of v_theta range but be > -pi and < pi.' );

Ntheta = Ntheta + 1;
v_theta = [ v_theta, pi ];
a_c_slice_radial_resampled = [ st_gen_resample.a_c_slice_radial_resampled; st_gen_resample.a_c_slice_radial_resampled(1,:) ];
a_dk_c_slice_radial_resampled = [ st_gen_resample.a_dk_c_slice_radial_resampled; st_gen_resample.a_dk_c_slice_radial_resampled(1,:) ];
a_dtheta_c_slice_radial_resampled = [ st_gen_resample.a_dtheta_c_slice_radial_resampled; st_gen_resample.a_dtheta_c_slice_radial_resampled(1,:) ];


% Sanitise the theta query vector to be in [-pi,pi)
[ v_theta_q ] = fn_ww__util__sanitise_theta_query_vec( v_theta_q );

v_c_slice_radial_resampled = reshape( a_c_slice_radial_resampled, 1, [] );
v_dk_c_slice_radial_resampled = reshape( a_dk_c_slice_radial_resampled, 1, [] );
v_dtheta_c_slice_radial_resampled = reshape( a_dtheta_c_slice_radial_resampled, 1, [] );                


if( v_theta(1) > v_theta(end))
    [~,~,v_theta_idxs] = histcounts(-v_theta_q,-v_theta);
else
    [~,~,v_theta_idxs] = histcounts(v_theta_q,v_theta);
end


assert( min( v_k ) > k_min && max( v_k ) < k_max, 'Interpolation ranges went iffy.' );

if( v_k(1) > v_k(end))
    [~,~,v_k_idxs] = histcounts(-v_k_q,-v_k);
else
    [~,~,v_k_idxs] = histcounts(v_k_q,v_k);
end     


% Optimisations (hopefully), CPU focused
l_v_k_idxs = length( v_k_idxs );
v_k_idxs_Ntheta = v_k_idxs * Ntheta;
v_k_idxs_m1_Ntheta = ( v_k_idxs - 1 ) * Ntheta;
v_k_idxs_Ntheta_pv_theta_idxs = v_k_idxs_Ntheta + v_theta_idxs;
v_k_idxs_m1_Ntheta_pv_theta_idxs = v_k_idxs_m1_Ntheta + v_theta_idxs;
v_k_idxs_Ntheta_pv_theta_idxs_p1 = v_k_idxs_Ntheta + v_theta_idxs + 1;
v_k_idxs_m1_Ntheta_pv_theta_idxs_p1 = v_k_idxs_m1_Ntheta + v_theta_idxs + 1;


% Do k radial interpolants, "only" cubic Hermite
v_c_radial_interp_A = zeros( 1, l_v_k_idxs );
v_c_radial_interp_B = zeros( 1, l_v_k_idxs );

v_k1 = v_k( v_k_idxs );
v_k2 = v_k( v_k_idxs+1 );
v_km = ( v_k1 + v_k2 ) / 2;
v_dk = v_k2 - v_km;
v_radial_s = ( v_k_q - v_km ) ./ v_dk;   % must be a row vector

a_radial_bf = [1/4*(v_radial_s+2).*(v_radial_s-1).^2
          1/4*(v_radial_s+1).*(v_radial_s-1).^2
          -1/4*(v_radial_s-2).*(v_radial_s+1).^2
          1/4*(v_radial_s-1).*(v_radial_s+1).^2];                          

v_c_radial_interp_A =  v_c_slice_radial_resampled( v_k_idxs_m1_Ntheta_pv_theta_idxs ) .* a_radial_bf(1,:) ...
                        + v_dk_c_slice_radial_resampled( v_k_idxs_m1_Ntheta_pv_theta_idxs ) .* a_radial_bf(2,:) .* v_dk ...
                        + v_c_slice_radial_resampled( v_k_idxs_Ntheta_pv_theta_idxs ) .* a_radial_bf(3,:) ...
                        + v_dk_c_slice_radial_resampled( v_k_idxs_Ntheta_pv_theta_idxs ) .* a_radial_bf(4,:) .* v_dk;
                    
v_c_radial_interp_B =  v_c_slice_radial_resampled( v_k_idxs_m1_Ntheta_pv_theta_idxs_p1 ) .* a_radial_bf(1,:) ...
                        + v_dk_c_slice_radial_resampled( v_k_idxs_m1_Ntheta_pv_theta_idxs_p1 ) .* a_radial_bf(2,:) .* v_dk ...
                        + v_c_slice_radial_resampled( v_k_idxs_Ntheta_pv_theta_idxs_p1 ) .* a_radial_bf(3,:) ...
                        + v_dk_c_slice_radial_resampled( v_k_idxs_Ntheta_pv_theta_idxs_p1 ) .* a_radial_bf(4,:) .* v_dk;

% Use c as weighting
c_A_weight = ( v_c_radial_interp_A - v_c_slice_radial_resampled( v_k_idxs_m1_Ntheta_pv_theta_idxs ) ) ...
                    ./ ( v_c_slice_radial_resampled( v_k_idxs_Ntheta_pv_theta_idxs ) - v_c_slice_radial_resampled( v_k_idxs_m1_Ntheta_pv_theta_idxs ) );

c_B_weight = ( v_c_radial_interp_B - v_c_slice_radial_resampled( v_k_idxs_m1_Ntheta_pv_theta_idxs_p1 ) ) ...
                    ./ ( v_c_slice_radial_resampled( v_k_idxs_Ntheta_pv_theta_idxs_p1 ) - v_c_slice_radial_resampled( v_k_idxs_m1_Ntheta_pv_theta_idxs_p1 ) );

dtheta_c_A = c_A_weight .* v_dtheta_c_slice_radial_resampled( v_k_idxs_m1_Ntheta_pv_theta_idxs ) ...
                    + ( 1 - c_A_weight ) .*  v_dtheta_c_slice_radial_resampled( v_k_idxs_Ntheta_pv_theta_idxs );

dtheta_c_B = c_B_weight .* v_dtheta_c_slice_radial_resampled( v_k_idxs_m1_Ntheta_pv_theta_idxs_p1 ) ...
                    + ( 1 - c_B_weight ) .*  v_dtheta_c_slice_radial_resampled( v_k_idxs_Ntheta_pv_theta_idxs_p1 );
 
               

% Do interpolant along angular direction
v_theta1 = v_theta( v_theta_idxs );
v_theta2 = v_theta( v_theta_idxs+1 );
v_thetam = ( v_theta1 + v_theta2 ) / 2;
v_dtheta = v_theta2 - v_thetam;        
v_angular_s = ( v_theta_q - v_thetam ) ./ v_dtheta;   % must be a row vector
      
a_angular_bf = [1/4*(v_angular_s+2).*(v_angular_s-1).^2
          1/4*(v_angular_s+1).*(v_angular_s-1).^2
          -1/4*(v_angular_s-2).*(v_angular_s+1).^2
          1/4*(v_angular_s-1).*(v_angular_s+1).^2];                      

v_c_PF_interp =  v_c_radial_interp_A .* a_angular_bf(1,:) ...
                    + dtheta_c_A .* a_angular_bf(2,:) .* v_dtheta ...
                    + v_c_radial_interp_B .* a_angular_bf(3,:) ...
                    + dtheta_c_B .* a_angular_bf(4,:) .* v_dtheta;  

       
        
end
function [ v_dkx_c_all, v_dky_c_all ] = fn_ww__util__grad_polar_to_cart( v_k_all, v_theta_all, v_dk_c_all, v_dtheta_c_all )
% fn_ww__util__grad_polar_to_cart: DEV

    v_dkx_c_all = v_dk_c_all .* cos( v_theta_all ) - ( 1 ./ v_k_all ) .* v_theta_all .* sin( v_theta_all );
    v_dky_c_all = v_dk_c_all .* sin( v_theta_all ) + ( 1 ./ v_k_all ) .* v_theta_all .* cos( v_theta_all );

end
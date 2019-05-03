function fn_ww__test__shear_profile__do_exp( Nk, k_min, k_max, Fr2, h, phy_U0 )
%fn_ww__test__shear_profile__do_exp: Test of powerlaw shear profile
%
%   fn_ww__test__shear_profile__do_exp( Nk, k_min, k_max, Fr2, h, phy_U0 )
%
% Test the exponential shear profile
%
% TAGS: CORE, WWERRINSHEAR
%
% See also
%   fn_ww__setup__shear_fn__nondim_exp(),
%   fn_ww__test__shear_profile__core(),
%   fn_ww__util__create_k_vec(),
%   fn_ww__setup__merge_parameters(),
%   fn_ww__test__shear_profile__do_powerlaw()



fprintf( 'Init...\n' );
[ st_p ] = fn_ww__setup__param_std__re_cl(  );
[ st_fn_shear_exp, st_p ] = fn_ww__setup__shear_fn__nondim_exp( st_p, Fr2, h, phy_U0 );
[ st_p_mp ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_mp', true ) );

[ v_k ] = fn_ww__util__create_k_vec( k_min, k_max, Nk, 3, 0 );


% Do tests
[ st_res ] = fn_ww__test__shear_profile__core( st_fn_shear_exp, v_k, k_min, k_max, st_p, st_p_mp );


% Simen test
fprintf( '... Simen... \n');
Fr = sqrt(Fr2);
v_c0_SE = 1./sqrt(v_k)/Fr; %Quiescent intrinsic velocity
zlen = 501;
a_kernel = zeros(Nk,zlen);
v_c_UDopp = v_c0_SE*0;
for i = 1:Nk
    zmin = max(-h,-3*pi/v_k(i));
    v_z_SE = zmin:-zmin/(zlen-1):0;
    v_U_SE = exp(v_z_SE);     %Velocity profile
    
	%This is 2k*cosh(2k(h+z))/sinh(2kh), re-written to manage big/small
	a_kernel = 2*v_k(i)*(exp(2*v_k(i)*v_z_SE)+exp(-2*v_k(i)*(2*h+v_z_SE)))/(1-exp(-4*v_k(i)*h));
    v_c_UDopp(i) = trapz(v_z_SE,v_U_SE.*a_kernel);
end
v_c_SE = v_c0_SE + v_c_UDopp;



% Note final parameters
st_res.Fr2 = Fr2;
st_res.Nk = Nk;
st_res.k_min = k_min;
st_res.k_max = k_max;
st_res.s_profile = 'exp';
st_res.phy_U0 = phy_U0;
st_res.h = h;
st_res.v_c_SE = v_c_SE;
st_res.st_p = st_p;

% Save it
save( 'dr_sim_exp', 'st_res' );


fprintf( '... done.\n' );








end
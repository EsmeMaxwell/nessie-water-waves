function fn_ww__test__shear_profile__do_columbia( Nk, k_min, k_max )
%fn_ww__test__shear_profile__do_columbia: Test of Columbia shear profile
%
%   fn_ww__test__shear_profile__do_columbia( Nk, k_min, k_max, Fr2, phy_h )
%
% Test the Colombia river shear profile.
%
% TAGS: CORE, WWERRINSHEAR
%
% See also
%   fn_ww__setup__shear_fn__nondim_powerlaw(),
%   fn_ww__test__shear_profile__core(),
%   fn_ww__util__create_k_vec(),
%   fn_ww__setup__merge_parameters(),
%   fn_ww__test__shear_profile__do_exp()




fprintf( 'Init...\n' );
[ st_p ] = fn_ww__setup__param_std__re_cl(  );
[ st_fn_shear_poly, st_p ] = fn_ww__setup__shear_fn__nondim_columbia_poly( st_p );
[ st_r_shear ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear_poly, st_p );
[ st_p_mp ] = fn_ww__setup__merge_parameters( st_p, struct( 'bp_mp', true ) );

[ v_k ] = fn_ww__util__create_k_vec( k_min, k_max, Nk, 3, 0 );


% Do tests
[ st_res ] = fn_ww__test__shear_profile__core( st_r_shear, v_k, k_min, k_max, st_p, st_p_mp );

% Note final parameters
st_res.Fr2 = st_p.Fr2;
st_res.Nk = Nk;
st_res.k_min = k_min;
st_res.k_max = k_max;
st_res.s_profile = 'columbia';
st_res.phy_U0 = st_p.phy_U0;
st_res.phy_h = st_p.phy_h;
st_res.st_p = st_p;


% Save it
save( 'dr_sim_columbia', 'st_res' );

fprintf( '... done.\n' );




end
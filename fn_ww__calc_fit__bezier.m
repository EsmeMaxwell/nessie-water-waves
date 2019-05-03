function [ st_fit_bezier ] = fn_ww__calc_fit__bezier( v_zm, v_zs, v_zs_candidate_U, h, hs, Nextrpts )
%fn_ww__calc_fit__bezier: Calc Bezier fit of data


% Check if surface point defined, if not then we must extrapolate
if ( 0 == hs )
    [ v_y, v_dy, v_ddy ] = fn_ww__dev__de_casteljau( v_zs_candidate_U, v_zm, -h, 0 );
else
    [ v_zs_candidate_U_extr ] = fn_ww__calc_fit__prep_lin_surf_extrapolate( v_zs, v_zs_candidate_U, Nextrpts );
    [ v_y, v_dy, v_ddy ] = fn_ww__dev__de_casteljau( v_zs_candidate_U_extr, v_zm, -h, 0 );
end

% Why do we need this?
v_y = flipud( v_y );
v_dy = flipud( v_dy );
v_ddy = flipud( v_ddy );

st_fit_bezier = struct;
st_fit_bezier.v_y = v_y;
st_fit_bezier.v_dy = v_dy;
st_fit_bezier.v_ddy = v_ddy;

% plot( v_zm, v_y )
% hold on;
% plot( v_zs, v_zs_candidate_U, 'r*' );
% hold off;
% error()


end
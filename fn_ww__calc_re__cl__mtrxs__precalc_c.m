function [ a_A0_precomp, a_A2, a_A0_FS, a_A1_FS, a_A2_FS ] = fn_ww__calc_re__cl__mtrxs__precalc_c( st_Dn, a_U, a_dU, a_ddU, a_I, st_p )
%fn_ww__calc_re__cl__mtrxs__precalc_c: Calc CL precalc matrices
% 
%   [ a_A0_precomp, a_A2, a_A0_FS, a_A1_FS, a_A2_FS ] = fn_ww__calc_re__cl__mtrxs__precalc_c( st_Dn, a_U, a_dU, a_ddU, a_I, st_p )
% 
% Calculates the precalc matrices for collocation. These only need be done
% once and can be reused (at least for the reduced problem).
% 
% INPUT
% 
%   st_Dn : Differentiation matrices
%   a_U, a_dU, a_ddU : Shear profile, diagonal matrices
%   a_I : Identity
%   st_p : Parameters
% 
% See also
%   fn_ww__calc_re__cl__mtrxs__core_c(),
%   fn_ww__calc_re__cl__red_c()


if ( st_p.bp_phy_calc )
    % Precomputations using diff matrices on phy domain
    a_A0_precomp = ( a_U * st_Dn.a_D2p - a_ddU  );
    a_A2 = 0 * a_I;

    % FS Precomputations using diff matrices on phy domain and g
    a_A0_FS = ( a_U^2 * st_Dn.a_Dp - a_dU * a_U - st_p.phy_g * a_I );
    a_A1_FS = ( -2 * a_U * st_Dn.a_Dp + a_dU );
    a_A2_FS = st_Dn.a_Dp; 
else
    % Precomputations using diff matrices on nondim domain
    a_A0_precomp = ( a_U * st_Dn.a_D2m - a_ddU  );
    a_A2 = 0 * a_I;

    % FS Precomputations using diff matrices on nondim domain and use the
    % nondim Froude
    a_A0_FS = ( a_U^2 * st_Dn.a_Dm - a_dU * a_U - ( 1 / st_p.Fr2 ) * a_I );
    a_A1_FS = ( -2 * a_U * st_Dn.a_Dm + a_dU );
    a_A2_FS = st_Dn.a_Dm; 
end

end
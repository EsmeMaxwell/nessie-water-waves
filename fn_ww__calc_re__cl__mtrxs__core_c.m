function [ a_A0, a_A1, a_A2 ] = fn_ww__calc_re__cl__mtrxs__core_c( st_Dn, a_U, a_A0_precomp, a_A2, a_I, a_A0_FS, a_A1_FS, a_A2_FS, k, st_p )
%fn_ww__calc_re__cl__mtrxs__core_c: Calc CL matrices A2,A1,A0 used in loop
% 
%   [ a_A0, a_A1, a_A2 ] = fn_ww__calc_re__cl__mtrxs__core_c( st_Dn, a_U, a_A0_precomp, a_A2, a_I, a_A0_FS, a_A1_FS, a_A2_FS, k, st_p )
% 
% Calculates the matrices used in the loop of the collocation calculation.
% Can handle both nondimensional and dimensional.
% 
% INPUT
%   st_Dn : Differentiation matrices
%   a_U : Shear profile, diagonal matrix
%   a_A0_precomp, a_A2, a_I, a_A0_FS, a_A1_FS, a_A2_FS : Precalculated
%   matrices
%   k : Wave number
%   st_p : Parameters
% 
% OUTPUT
%
%   a_A0, a_A1, a_A2 : Matrices used for CL calc
% 
% See also
%   fn_ww__calc_re__cl__mtrxs__precalc_c(),
%   fn_ww__calc_re__cl__red_c()




if ( st_p.bp_phy_calc )
    % Main equations using diff matrices on phys domain
    a_A0 = ( a_A0_precomp - k^2 * a_U );
    a_A1 = -st_Dn.a_D2p + k^2 * a_I;

    a_A0(1,:) = a_A0_FS(1,:);
    a_A1(1,:) = a_A1_FS(1,:);
    a_A2(1,:) = a_A2_FS(1,:);

    a_A0(end,:) = 0 * a_A0(end,:);
    a_A1(end,:) = 0 * a_A1(end,:);
    a_A2(end,:) = 0 * a_A2(end,:);

else
    % Main equations using diff matrices on nondim domain
    a_A0 = ( a_A0_precomp - k^2 * a_U );
    a_A1 = -st_Dn.a_D2m + k^2 * a_I;

    a_A0(1,:) = a_A0_FS(1,:);
    a_A1(1,:) = a_A1_FS(1,:);
    a_A2(1,:) = a_A2_FS(1,:);

    a_A0(end,:) = 0 * a_A0(end,:);
    a_A1(end,:) = 0 * a_A1(end,:);
    a_A2(end,:) = 0 * a_A2(end,:);

end


end
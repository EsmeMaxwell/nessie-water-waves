function [ a_A2, a_A1, a_A0, alpha ] = fn_ww__util__fld_scaling( a_A2, a_A1, a_A0 )
%fn_ww__util__fld_scaling: Util FLD scaling
%
%   [ a_A2, a_A1, a_A0, alpha ] = fn_ww__util__fld_scaling( a_A2, a_A1, a_A0 )
%
% Perform scaling of matrices in quadtratic eigenvalue problem. Imroves,
% albeit marginally in our case, the eigenvalue resutls but may cause
% issues when calculating backwards error and similar.
%
% See fn_ww__setup__param_std__re_cl() and parameter bp_fld.
%
% SEE : Hung-Yuan Fan, Wen-Wei Lin, and Paul Van Dooren, "Normwise Scaling
% of Second Order Polynomial Matrices", SIAM J. Matrix Anal. Appl., 26(1),
% 252–256
%
% URL : https://epubs.siam.org/doi/abs/10.1137/S0895479803434914?journalCode=sjmael
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__calc_re__cl__red_c()


gamma0 = norm( a_A0, 2 );
gamma1 = norm( a_A1, 2 );
gamma2 = norm( a_A2, 2 );
alpha = sqrt( gamma0 / gamma2 );
beta = 2 / ( gamma0 + gamma1 * alpha );

a_A2 = alpha^2 * beta * a_A2;
a_A1 = alpha * beta * a_A1;
a_A0 = beta * a_A0;

end
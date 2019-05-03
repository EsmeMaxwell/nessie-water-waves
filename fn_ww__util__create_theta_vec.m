function [ v_theta ] = fn_ww__util__create_theta_vec( Ntheta )
%fn_ww__util__create_theta_vec: Util creates uniformally spaces theta vector
% 
%   [ v_theta ] = fn_ww__util__create_theta_vec( Ntheta )
% 
% Create a theta vector properly, with uniformally spaces angles.
% 
% INPUT
%   Ntheta : Number of elements required
% 
% OUTPUT
%   v_theta : Vector of theta angles
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__util__create_k_vec(),
%   fn_ww__calc_re__cl__gen_sl_pol_ang_c()

v_theta = ( 2 * pi * ( 0:1:(Ntheta-1) ) / Ntheta ) - pi;

end
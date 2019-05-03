function [ v_theta_q ] = fn_ww__util__sanitise_theta_query_vec( v_theta_q )
%fn_ww__util__sanitise_theta_query_vec: Util sanitises a theta vector to [0,2\pi)
% 
%   [ v_theta_q ] = fn_ww__util__sanitise_theta_query_vec( v_theta_q )
% 
% Remaps points into interval [0,2\pi) assuming 2\pi periodicity.
% 
% INPUT
%   v_theta : Vector of theta angles
%
% OUTPUT 
%   v_theta : Vector of theta angles in [0,2\pi)
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
% 

v_theta_q = mod( v_theta_q + pi, 2*pi ) - pi;

end
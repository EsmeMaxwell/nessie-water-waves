function [tc_kc,tc_EL] = fn_ww__ext__dim__calc_KC_EL(I,Ik,z,b_no_EL)
%fn_ww__ext__dim__calc_KC_EL: EXTERNAL Yan Li DIM calc KC and EL approximations
% 
%   [tc_kc,tc_EL] = fn_ww__ext__dim__calc_KC_EL(I,Ik,z,b_no_EL)
% 
% Calculate KC and EL approximations; this is a modified version of Yan
% Li's DIM code.
%
% TAGS: EXT
%
% See also
%   fn_ww__ext__dim__c()


%%  inputs
g       = I.g;
h       = I.h;
kv      = Ik.k;
k_z     = Ik.k_i;
theta   = Ik.tht;
kx_z    = k_z.*cos(theta);
ky_z    = k_z.*sin(theta);
z_k     = z.z_k;
dz      = z.dz; %h/I.N
   
%%       
th      = tanh(kv.*h);
c0      = sqrt(g.*th./kv);
   
%% E&L


%% 
sh2sh    = (exp(2*k_z.*z_k)-exp(-4*k_z*h-2*k_z.*z_k))./(1-exp(-4.*k_z*h));
itg_kc   = (kx_z.*I.dUkx_i+ky_z.*I.dUky_i)./k_z.*sh2sh;

tc_kc    = c0-fn_ww__ext__dim__simpsons(itg_kc,dz);

if ( b_no_EL ~= 1 )
    dltL     = fn_ww__ext__dim__simpsons(itg_kc,dz)./c0;
    tc_EL    = c0.*(sqrt(1+dltL.^2)-dltL);
else
    tc_EL = 0 * tc_kc;
end



end


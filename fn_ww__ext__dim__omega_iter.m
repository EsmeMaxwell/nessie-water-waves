function [tc,wz,R_err0,itr_cnt] = fn_ww__ext__dim__omega_iter(I,Ik,tc0,z,tol,M)
%fn_ww__ext__dim__omega_iter: EXTERNAL Yan Li DIM calc iteratative scheme
% 
%   [tc,wz,R_err0,itr_cnt] = fn_ww__ext__dim__omega_iter(I,Ik,tc0,z,tol,M)
% 
% Calculate iteratative scheme for DIM; this is Yan Li's DIM code.
%
% TAGS: EXT
%
% ORIGINAL HEADERS
%
% either tol or M (the iterative number) controls the
% convergence
%
% return: tc: the relative phase velocity; wz:the vertical velocity;
% R_err0: the maximum relative error
%
%
% See also
%   fn_ww__ext__dim__c()



%% size(k)*N
g       = I.g;
h       = I.h;
k_z     = Ik.k_i;
k       = Ik.k;
kx      = Ik.kx;
ky      = Ik.ky;
kx_z    = Ik.kx_i;
ky_z    = Ik.ky_i;
z_k     = z.z_k;
dz      = z.dz;  % now a vector of the same size of k; PM is the k-dependent z spacing


% - % -% - % - % - % - %
%     inputs above     %
% - % -% - % - % - % - % 

%%
thdk    = tanh(k*h)./k;
kdu0    = kx*I.dU0(1)+ky*I.dU0(2);

ku0_z   = (kx_z*I.U0(1)+ky_z*I.U0(2));
kuz_z   = (kx_z.*I.Ukx_i+ky_z.*I.Uky_i);
kDtu0_z =  kuz_z-ku0_z;                 % PM: k.\Delta U
kdduz_z = (kx_z.*I.ddUkx_i+ky_z.*I.ddUky_i); 
kduz    = (kx_z.*I.dUkx_i+ky_z.*I.dUky_i);
shdch   = (exp(k_z.*z_k)-exp(-k_z.*z_k-2*k_z*h))./(1+exp(-2*k_z*h));

%% obtain some parameters from the Rayleight Equation as a function of k and c
tc           = tc0.tc;
tc_z         = tc0.tc_i;
ku_minus_kcz = kDtu0_z-k_z.*tc_z;           % PM: k.\Delta U - kc  (in a z vector)
kddu_dv_kc   = kdduz_z./ku_minus_kcz;       % PM: ( k.U'' ) / ( k.\Delta U - kc )
kddu_dv_kc(isnan(kddu_dv_kc))=0;  % avoid a critical layer

%% estimate the initial error
kbcz      = kx_z.^2+ky_z.^2+ kddu_dv_kc;        % PM: Rayleigh eqn when rearrangmed for w'' (and no post w mult)
[wz]      = fn_ww__ext__dim__get_w(I.N,dz,kbcz);  % obtain w when an inital guess of c is used
ig_itg    = kdduz_z.*wz.*shdch./k_z./ku_minus_kcz;
dig_itg   = kdduz_z.*wz.*shdch./ku_minus_kcz.^2;
ig_itg(isnan(ig_itg))  = 0;
dig_itg(isnan(dig_itg))= 0;

d_ig      = fn_ww__ext__dim__simpsons(dig_itg,dz);
ig        = fn_ww__ext__dim__simpsons(ig_itg,dz);
delta_c   = (1+ig).*tc.^2+tc.*thdk.*kdu0./k-g.*thdk;
d_delta_c = 2*(1+ig).*tc+thdk.*kdu0./k+d_ig.*tc.^2;
converg   = abs(delta_c./tc./d_delta_c);  % the relative error defined
R_err0    = converg;
max_convg = max(converg);
j         = M;  %iterative number

itr_cnt = 0;

%% start a loop to calculate an exact solution
while max_convg>tol 
    itr_cnt = itr_cnt + 1;
    j            = j-1;
    epw          = 1; 
    tc           = tc-epw.*delta_c./d_delta_c;    % PM: Newton iteration?, x_{n+1} = x_n - f_n / f_n'

    ku_minus_kcz = kDtu0_z-k_z.*tc;
    kddu_dv_kc   = kdduz_z./ku_minus_kcz;
    
    kddu_dv_kc(isnan(kddu_dv_kc))=0;  % avoid a critical layer   
    
    kbcz         = kx_z.^2+ky_z.^2+ kddu_dv_kc;
    [wz]         = fn_ww__ext__dim__get_w(I.N,dz,kbcz);
    ig_itg       = kdduz_z.*wz.*shdch./k_z./ku_minus_kcz;
    dig_itg      = kdduz_z.*wz.*shdch./ku_minus_kcz.^2;
    
    ig_itg( isnan(ig_itg) )  = 0;
    dig_itg( isnan(dig_itg) )= 0;
    
    d_ig         = fn_ww__ext__dim__simpsons( dig_itg, dz );
    ig           = fn_ww__ext__dim__simpsons( ig_itg, dz );
    delta_c      = (1+ig).*tc.^2+tc.*thdk.*kdu0./k-g.*thdk;
    d_delta_c    = 2*(1+ig).*tc+thdk.*kdu0./k+d_ig.*tc.^2;
    converg      = abs(delta_c./tc./d_delta_c);
    max_convg    = max(converg);
    if j<1
        break
    end
end


end

function [tc_app]= c_app(g,h,kv,kduz,zv,dz)
[k,z]   = meshgrid(kv,zv);
c0      = sqrt(g./kv.*tanh(kv.*h));
nkc0    = sqrt(g.*k.*tanh(k.*h));

%%
sh2sh   = (exp(2.*k.*z)-exp(-2.*k.*z-4.*k.*h))./(1-exp(-4.*k.*h));
dltL    = fn_ww__ext__dim__simpsons(kduz.*sh2sh./nkc0,dz);
tc_app  = c0.*(sqrt(1+dltL.^2)-dltL);
% tc_kc   = c0.*(1-dltL);

end
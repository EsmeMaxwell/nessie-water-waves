function [ tc_kc, tc_EL, tc_dimM, itr_cnt_short, itr_cnt_long ] = fn_ww__ext__dim__c( v_k, theta, st_fn, tol, N, max_itr, h, g, b_no_EL )
%fn_ww__ext__dim__c: EXTERNAL Yan Li DIM calc dispersion relation
% 
%   [ tc_kc, tc_EL, tc_dimM, itr_cnt_short, itr_cnt_long ] = fn_ww__ext__dim__c( v_k, theta, st_fn, tol, N, max_itr, h, g, b_no_EL )
% 
% Calculate dispersion relation; this is a modified version of Yan Li's DIM
% code.
% 
% TAGS: EXT
% 

%N = 2000;  % the total number of points in the vertical direction

%g = 1;
%h = 1;
%theta = 0;   
%tol       = 10e-10;  % may be used to control accuracy/convergence



% define short and long waves
nk        = length(v_k);            % number of k_vec
j = nk;
kh_cut    = ( 3.5 + 2*log(N/4) );  % used to cut z = -kh_cut/k to z = -h for short waves/deep waters
% kh_cut    = 90;
for i     = 1:nk    
    if v_k(i)> kh_cut
        j = i-1;
        break
    end
end

v_k_L     = v_k(1:j);        %long waves range




if ( j < nk )
        
    v_k_S     = v_k(j+1:end);    %short waves range


    % now short waves
    dz        = kh_cut./abs(v_k_S)./N; % for the DIM; PM is k-dependent z spacing
    z_cut     = -kh_cut./abs(v_k_S);
    dm_nz     = zeros(N+1,nk-j);  % z meshes; PM Nz rows, Nk columns


    for i=1:(nk-j)
        dm_nz(:,i)= 0:-dz(i):-kh_cut/abs(v_k_S(i));
    end


    % to obtain input information below
    Ik.k        = v_k_S;
    Ik.kx       = v_k_S.*cos(theta);  
    Ik.ky       = v_k_S.*sin(theta);
    [Ik.k_i,~]  = meshgrid(v_k_S,dm_nz(:,1));
    [Ik.kx_i,~] = meshgrid(Ik.kx,dm_nz(:,1)); 
    Ik.ky_i     = sign(Ik.ky).*sqrt(Ik.k_i.^2-Ik.kx_i.^2); 
    %[I,z]       = shearprofile(px,py,dm_nz);
    [I,z]       = fn_ww__ext__dim__setup_shear_fn( dm_nz, st_fn, h, g  );
    z.dz        = dz; 
    z.dm_nz     = dm_nz;

    I.N       = N; % intervals, N+1 nodes      
    I.g       = g;
    I.h       = h;
    Ik.tht    = theta;

    % to obtain KC and EL: tc = c - k_vec.*U_vec_0./|k_vec|; 
    [tc_kcS,tc_ELS]    = fn_ww__ext__dim__calc_KC_EL(I,Ik,z,b_no_EL); % obtan results from KC and EL

    tc0.tc             = tc_kcS;
    [tc0.tc_i,~]       = meshgrid(tc_kcS,dm_nz(:,1));  
    [tc_dimM_S,~,~,itr_cnt_short]    = fn_ww__ext__dim__omega_iter(I,Ik,tc0,z,tol,max_itr); % obtain results from DIM, in the input 1 denotes pnly one step of iteration

else
    
    % We don't bother calculating short waves here.
    tc_kcS = [];
    tc_ELS = [];
    tc_dimM_S = [];
    itr_cnt_short = 0;
    
end





% Long waves
dz        = h/N; % for the DIM 
z_v       = 0:-dz:-h;
[~,dm_nz] = meshgrid(v_k_L,z_v);



Ik.k        = v_k_L;
Ik.kx       = v_k_L.*cos(theta);  
Ik.ky       = v_k_L.*sin(theta);
[Ik.k_i,~]  = meshgrid(v_k_L,dm_nz(:,1));
[Ik.kx_i,~] = meshgrid(Ik.kx,dm_nz(:,1)); 
Ik.ky_i     = sign(Ik.ky).*sqrt(Ik.k_i.^2-Ik.kx_i.^2); 

%[I,z]       = shearprofile(px,py,dm_nz);
[I,z]       = fn_ww__ext__dim__setup_shear_fn( dm_nz, st_fn, h, g  );

    I.g       = g;
    I.h       = h;
    Ik.tht    = theta;

I.N         = N;
z.dm_nz     = dm_nz;
z.dz        = dz; % h/N
[tc_kcL,tc_ELL]    = fn_ww__ext__dim__calc_KC_EL(I,Ik,z,b_no_EL);

tc0.tc             = tc_kcL;
[tc0.tc_i,~]       = meshgrid(tc_kcL,dm_nz(:,1));
[tc_dimM_L,~,~,itr_cnt_long]    = fn_ww__ext__dim__omega_iter(I,Ik,tc0,z,tol,max_itr);



tc_kc       = [tc_kcL,tc_kcS];
tc_EL       = [tc_ELL,tc_ELS];
tc_dimM     = [tc_dimM_L,tc_dimM_S];

% % TEMP change
tc_kc = tc_kc + I.U0(1)*cos(theta) + I.U0(2)*sin(theta);
tc_EL = tc_EL + I.U0(1)*cos(theta) + I.U0(2)*sin(theta);
tc_dimM = tc_dimM + I.U0(1)*cos(theta) + I.U0(2)*sin(theta);



end
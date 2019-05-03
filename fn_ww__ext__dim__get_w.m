function [nwz] = fn_ww__ext__dim__get_w(n,dz,kbcz)
%fn_ww__ext__dim__calc_KC_EL: EXTERNAL Yan Li DIM calc w vector from Dirichlet conditions
% 
%   [nwz] = fn_ww__ext__dim__get_w(n,dz,kbcz)
% 
% Calculate w vector using Dirichlet BCs; this is Yan Li's code.
%
% TAGS: EXT
%
% See also
%   fn_ww__ext__dim__c()




%%      can develop to nk*nz matrix, so that the speed can be highly improved
nwz      = 0.*kbcz;
nk       = length(kbcz(1,:));
% nwz      = zeros(1,n+1);
nwz(1,:) = 1;

%%
%----------------------------------------%
%---- obtain the coefficient matri   ----%
%----------------------------------------%
cfM_11   = dz.^2.*kbcz(2,:)+2;  %w2
cfM_12   = -1;

%%
rht      = zeros(n-1,nk);  
y        = rht;
gm       = y;

rht(1,:) = nwz(1,:);

%----------------------------------------%
%--- start the iteration for the matrix--%
%----------------------------------------%
%% b_i diagonal, a_i, lower half, c_i upper half
y(1,:)  = rht(1,:)./cfM_11;
gm(1,:) = cfM_12./cfM_11;
for  i  = 2:(n-2)
    %% from governing equation
        cfM_ii   = dz.^2.*kbcz(i+1,:)+2;  
%         cfM_iip1 = -1;   
%         cfM_iim1 = -1;
%         aph_i    = cfM_iim1;
%         beta_i   = cfM_ii-aph_i.*gm(i-1,:);
%         gm(i,:)  = cfM_iip1./beta_i;
        beta_i   = cfM_ii+gm(i-1,:);
        gm(i,:)  = -1./beta_i;
        y(i,:)   = (rht(i,:)+y(i-1,:))./beta_i;        
end

%% w_n
cfM_ii    = dz.^2.*kbcz(n,:)+2;  
% cfM_iim1  = -1;
% i         = n-1;
% aph_i     = cfM_iim1;
% beta_i    = cfM_ii-aph_i.*gm(i-1,:); 
% y(i,:)      = (rht(i,:)-aph_i.*y(i-1,:))./beta_i;  
i           = n-1;
beta_i      = cfM_ii+gm(i-1,:); 
y(i,:)      = (rht(i,:)+y(i-1,:))./beta_i;   
nwz(n,:)    = y(n-1,:);
for i = (n-2):-1:1
    nwz(i+1,:) = y(i,:)-gm(i,:).*nwz(i+2,:);
end
% wz(2:1:n) = rts(1:(n-1));


end
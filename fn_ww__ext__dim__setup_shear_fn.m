function [I,z] = fn_ww__ext__dim__setup_shear_fn( dm_nz, st_fn, h, g )
%fn_ww__ext__dim__setup_shear_fn: EXTERNAL Yan Li DIM setup shear profile
% 
%   [I,z] = fn_ww__ext__dim__setup_shear_fn( dm_nz, st_fn, h, g )
% 
% Setup shear profile; this is a modified version of Yan Li's DIM code.
%
% TAGS: EXT
%
% See also
%   fn_ww__ext__dim__create_dim_shear(),
%   fn_ww__ext__dim__c()




%% define the shear profile   
   %  pxinput is the coefficients vector of Ux(z)
   %  pyinput is the coefficients vector of Uy(z)
%         p1       = pxinput(1);
%         p2       = pxinput(2);
%         p3       = pxinput(3);
%         p4       = pxinput(4);
%         p5       = pxinput(5);
%         p6       = pxinput(6);
% %         p7       = pxinput(7);
%         p1y      = pyinput(1);
%         p2y      = pyinput(2);
%         p3y      = pyinput(3);
%         p4y      = pyinput(4);
%         p5y      = pyinput(5);
%         p6y      = pyinput(6);
% %         p7y      = pyinput(7);
%         
%         uxf      = @(z) polyval(pxinput,z);
%         pcofdx   = [6*p1 5*p2 4*p3 3*p4 2*p5 p6];
%         duxf     = @(z) polyval(pcofdx,z);       
%         pcofddx  = [30*p1 20*p2 12*p3 6*p4 2*p5];
%         dduxf    = @(z) polyval(pcofddx,z);
%         
%         uyf      = @(z) polyval(pyinput,z);       
%         pcofdy   = [6*p1y 5*p2y 4*p3y 3*p4y 2*p5y p6y];
%         duyf     = @(z) polyval(pcofdy,z);        
%         pcofddy  = [30*p1y 20*p2y 12*p3y 6*p4y 2*p5y];
%         dduyf    = @(z) polyval(pcofddy,z);
%         

% % for h=2 tests
% csp_a = 0.3;
% csp_b = 0.2;

% for h=10 tests

% csp_a = 0.01;
% csp_b = 0.15;
% 
%         % A simple(ish) profile
%         uxf = @(z) csp_b * ( cos( csp_a * (z).^2 ) );
%         duxf = @(z) csp_b * ( -2 * csp_a * (z) .* sin( csp_a * (z).^2 ) );
%         dduxf = @(z) csp_b * ( -4 * csp_a^2 * (z).^2 .* cos( csp_a * (z).^2 ) -2 * csp_a * sin( csp_a * (z).^2 ) );


% if ( isstruct( st_U_pp ) ) 
% 
%     uxf = @(z) ppval( st_U_pp, z );   
%     duxf = @(z) ppval( st_dU_pp, z );   
%     dduxf = @(z) ppval( st_ddU_pp, z );   
%         
% else
% 
%     uxf = @(z) polyval( st_U_pp, z );   
%     duxf = @(z) polyval( st_dU_pp, z );   
%     dduxf = @(z) polyval( st_ddU_pp, z );       
%     
% end
%         uyf = @(z) 0 * z;
%         duyf = @(z) 0 * z;
%         dduyf = @(z) 0 * z;
%         

        uxf = st_fn.fn_Ux;
        duxf = st_fn.fn_dUx;
        dduxf = st_fn.fn_ddUx;
        
        uyf = st_fn.fn_Uy;
        duyf = st_fn.fn_dUy;
        dduyf = st_fn.fn_ddUy;

        
        
        %% VALUES
%         I.U_i(1,:)   = uxf(dm_nz);
%         I.U_i(2,:)   = uyf(dm_nz);
%         I.dU_i(1,:)  = duxf(dm_nz);
%         I.dU_i(2,:)  = duyf(dm_nz);
%         I.ddU_i(1,:) = dduxf(dm_nz);
%         I.ddU_i(2,:) = dduyf(dm_nz);        
        I.U_i        = [];
        I.dU_i       = [];
        I.ddU_i      = [];    
        I.U0(1)      = uxf(0);
        I.U0(2)      = uyf(0);
        I.dU0(1)     = duxf(0);
        I.dU0(2)     = duyf(0);        
        
        %% place holders
        I.g       = g;
        I.h       = h;
        I.N       = []; % intervals, N+1 nodes
        z.n_dz    = []; %1*nz 
        z.dz      = []; % h/N
       
        
        %% values
        z.dm_nz  = []; % vector 1*nz;
        z.z_k    = dm_nz;
        
        %% size = nz*nk
        [I.Ukx_i]      = uxf(dm_nz);
        [I.Uky_i]      = uyf(dm_nz);
        [I.dUkx_i]     = duxf(dm_nz);
        [I.dUky_i]     = duyf(dm_nz);
        [I.ddUkx_i]    = dduxf(dm_nz);
        [I.ddUky_i]    = dduyf(dm_nz);
            
end
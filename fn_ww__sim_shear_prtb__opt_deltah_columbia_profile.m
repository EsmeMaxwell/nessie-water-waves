function fn_ww__sim_shear_prtb__opt_deltah_columbia_profile(  )
%fn_ww__sim_shear_prtb__full_columbia_profile: Sim manager to run err-in-shear for Columbia profile to find optimal delta h
%
%   fn_ww__sim_shear_prtb__opt_deltah_columbia_profile(  )
%
% Manager function to setup and execute the err-in-shear-profile simulation
% for the columbiaonential profile, this time to find optimal delta h.
% 
% See also
%   fn_ww__sim_shear_prtb__do_runs(),
%   fn_ww__setup__shear_fn__nondim_columbia()


Nitr = 1000;
Nz = 96;
Nk = 150;
phy_h = 20;
phy_U0 = 1;
v_percentile = 5:5:95;



% Setup parameter space
[ st_p ] = fn_ww__setup__param_std__re_cl(  );

% Setup profile and update parameter set    
[ st_fn_shear_poly, st_p ] = fn_ww__setup__shear_fn__nondim_columbia_poly( st_p );

% Create log distributed k vector
[ v_k ] = fn_ww__util__create_k_vec( 0.3, 300, Nk, 1, 0 );

% Simulation parameter setup
v_Fr2 = [ st_p.Fr2 ];
v_hs_prop = [ 0 ];
%v_delta_h_prop = 0.0025:0.0025:0.08;
v_delta_U = [ 0.1 ];

Ndh = 20;
v_delta_h_prop_initial = fn_ww__util__create_k_vec( (1/250), 0.1, Ndh, 3, 0 );
v_required = union( 10:2:30, [ 10 25 30 35 40 45 50 55 60 65 75 100 125 250 500 ] );
v_delta_h_prop = 1 ./ unique( union( round( 1 ./ v_delta_h_prop_initial ), v_required ) )



NFr2 = numel( v_Fr2 );
Nhs = numel( v_hs_prop );
Ndeltah = numel( v_delta_h_prop );
NdeltaU = numel( v_delta_U );

Nspace = NFr2 * Nhs * Ndeltah * NdeltaU;
fprintf( '\n\nAbout to do %d test runs at %d iterations per run...\n\n', Nspace, Nitr );




st_param_save = struct;
st_param_save.v_Fr2 = v_Fr2;
st_param_save.v_hs_prop = v_hs_prop;
st_param_save.v_delta_h_prop = v_delta_h_prop;
st_param_save.v_delta_U = v_delta_U;
st_param_save.Nitr = Nitr;
st_param_save.Nz = Nz;
st_param_save.Nk = Nk;

ca_param_space = cell( NFr2, Nhs, Ndeltah, NdeltaU );



% Update parameter vectors correctly (usually this in the Fr2 loop)
v_hs = st_p.h * v_hs_prop;
v_delta_h = st_p.h * v_delta_h_prop;

% Oh, yeah. Coding like it's ninenteen ninety nine...
i_Fr2 = 1;
i_counter = 0;
    
% Now do the standard simulation...
for lp_hs=1:Nhs
    for lp_dh=1:Ndeltah
        for lp_dU=1:NdeltaU

            i_counter = i_counter + 1;
            fprintf( 'Doing run %d of %d (%d,%d,%d,%d): Fr2=%0.2f; hs=%0.2f; dh=%0.2f; dU=%0.2f.\n', i_counter, Nspace, i_Fr2, lp_hs, lp_dh, lp_dU, v_Fr2(i_Fr2), v_hs(lp_hs), v_delta_h(lp_dh), v_delta_U(lp_dU) );               

            % Setup measurement vectors
            h_usable = st_p.h - v_hs(lp_hs);
            h_pts = floor( h_usable / v_delta_h(lp_dh) );
            v_zs = linspace( -v_hs(lp_hs), -st_p.h, h_pts ).';
            v_zs_err_width = v_delta_U(lp_dU) * ones( size( v_zs ) );
            fprintf( 'Delta h = %f; hs = %f, num pts = %d\n', v_delta_h(lp_dh), v_hs(lp_hs), h_pts );

            % Create suitable filename
            s_filename = sprintf( 'data_shear_prtb_opt_deltah_columbia/shear_prtb_columbia_opt_deltah__%d_%d_%d_%d', i_Fr2, lp_hs, lp_dh, lp_dU );
            s_backup_filename = sprintf( 'data_shear_prtb_opt_deltah_columbia/shear_prtb_columbia_opt_deltah__data__%d_%d_%d_%d', i_Fr2, lp_hs, lp_dh, lp_dU );

            % Save parameter space details
            st_param = struct;
            st_param.s_filename = s_filename;
            st_param.v_hs_prop = v_hs_prop;
            st_param.v_delta_h_prop = v_delta_h_prop;
            st_param.Fr2 = v_Fr2(i_Fr2);
            st_param.hs = v_hs(lp_hs);
            st_param.delta_h = v_delta_h(lp_dh);
            st_param.delta_U = v_delta_U(lp_dU);
            st_param.v_zs = v_zs;
            st_param.v_zs_err_width = v_zs_err_width;
            st_param.h_pts = h_pts;
            st_param.st_p = st_p;
            ca_param_space{i_Fr2,lp_hs,lp_dh,lp_dU} = st_param;

            % Run the test
            tic
            [ st_err_run ] = fn_ww__sim_shear_prtb__do_runs( st_fn_shear_poly, v_zs, v_zs_err_width, v_k, v_percentile, Nitr, Nz, st_p, s_filename, s_backup_filename );
            toc;

            if ( ~st_err_run.b_ok )
                warning( 'Error with run' );
            end

            % Save the parameter space just incase we have to quit half
            % way through
            st_param_save.ca_param_space = ca_param_space;               
            save( 'data_shear_prtb_opt_deltah_columbia/shear_prtb_columbia_opt_deltah__param_space', 'st_param_save' );

        end
    end
end


end
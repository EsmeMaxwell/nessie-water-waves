function [ v_c_p, a_c_w_p ] = fn_ww__calc_re__adptv__red_c( st_Dn, st_r_shear, st_adptv, st_proc, st_p )
%fn_ww__calc_re__cl_adptv__red_c: Calc CL-c with adaptive depth and remapping of eigenvector for reduced profile
% 
%   [ v_c_p, a_c_w_p ] = fn_ww__calc_re__adptv__red_c( st_Dn, v_k, st_r_shear, st_proc, st_p )
% 
% 
% Must precalc data using fn_ww__calc_re__adptv_interval_refine()
% 
% Note, can probably make this code more efficient by paying some attention
% to MATLAB's copy-on-write semantics.

v_k = st_adptv.v_k;

Nk = numel( v_k );
pf_extra = 1e-5;

% Whether to use CL or PF, and whether to calc eigenvector
ip_method = st_proc.ip_method;
if ( nargout > 1 )
    bp_calc_eigvec = true;    
    v_zm_dst = st_proc.v_zm_dst;
    Nzdst = numel( v_zm_dst );
else
    bp_calc_eigvec = false;
end


% % Calculate intervals, optimal h, and overlap data
% [ st_adtpv ] = fn_ww__calc_re__adptv_interval_refine( v_k, st_p );

% For each interval, calculate that portion of the dispersion relation
% curve
ca_dr = cell( 1, st_adptv.Nintervals );
for lp_i=1:st_adptv.Nintervals

    %fprintf( '\nInterval %d: [ %f, %f ]...\n', lp_i, v_intervals(lp_i,1), v_intervals(lp_i,2) );
    
    % Local k interval
    v_k_loc = st_adptv.ca_v_k{lp_i}.v_k;
    
    % Adjust depth
    st_p_loc = st_p;
    st_p_loc.h = st_adptv.v_optimal_h(lp_i);
    st_p_loc.a = -st_adptv.v_optimal_h(lp_i);
    st_p_loc.bp_disp_update = false;    % TODO use the merge function!
    [ st_Dn_loc ] = fn_ww__setup__lin_map_Dn_to_mapped( st_Dn, st_p_loc );
        
    
    if ( 1 == ip_method )        
        
        % Calculate DR using CL
        if ( bp_calc_eigvec ) 
            [ v_c_p_loc, a_c_w_p_loc ] = fn_ww__calc_re__cl__red_c( st_Dn_loc, v_k_loc, st_r_shear, st_p_loc );
        else
            [ v_c_p_loc ] = fn_ww__calc_re__cl__red_c( st_Dn_loc, v_k_loc, st_r_shear, st_p_loc );
        end
        
    elseif ( 2 == ip_method )        
        
        % Calculate DR using PF
        k0 = 0.5 * ( v_k_loc(1) + v_k_loc(end ) );        

        [ st_dp ] = fn_ww__calc_re__clpf_dp__red_c( st_Dn_loc, v_k_loc(1)-pf_extra, v_k_loc(end)+pf_extra, st_r_shear, k0, st_proc.pf_tol, st_p_loc );
        
        if ( bp_calc_eigvec ) 
            [ a_y0_PF_interp ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp, v_k_loc, true );
            a_c_w_p_loc = a_y0_PF_interp( 1:end-1, : );
        else
            [ a_y0_PF_interp ] = fn_ww__calc_re__clpf_interp__sl_rad_c( st_dp, v_k_loc, false );            
        end
        
        v_c_p_loc = a_y0_PF_interp( end, : );
        
    else        
        error( 'Must specify valid computational method' );
    end
    
    % Remap eigenvector and savaeremapped c and w for this interval
	ca_dr{lp_i}.v_c_p = v_c_p_loc;
    if ( bp_calc_eigvec )     
        [ a_w_dst ] = fn_ww__interp__extend_short_eigenvector( st_Dn_loc, a_c_w_p_loc, st_p_loc, v_zm_dst, st_p );    
        ca_dr{lp_i}.a_w_p = a_w_dst;
    end
    
end




v_c_p = zeros( 1, Nk );
v_c_p( 1:numel(ca_dr{1}.v_c_p) ) = ca_dr{1}.v_c_p;

if ( bp_calc_eigvec )
    a_c_w_p = zeros( Nzdst, Nk );
    a_c_w_p( :, 1:size( ca_dr{1}.a_w_p, 2 ) ) = ca_dr{1}.a_w_p;
end
   
for lp_i=1:st_adptv.Noverlaps

    % Number of weight entries
    Nwrt = numel( st_adptv.ca_overlap{lp_i}.v_wt_r );
    
    % Do eigenvalue first (so code is less messy with conditional for
    % eigenvector)
        
    % Weight the end of the main interval using weights from overlap LHS
    v_c_p( st_adptv.ca_overlap{lp_i}.v_k_idxs ) = v_c_p( st_adptv.ca_overlap{lp_i}.v_k_idxs ) .* st_adptv.ca_overlap{lp_i}.v_wt_l;
    
    % Get next interval
    v_c_p_temp = ca_dr{lp_i+1}.v_c_p;
    
    % Weight 
    v_c_p_temp( 1:Nwrt ) = v_c_p_temp( 1:Nwrt ) .* st_adptv.ca_overlap{lp_i}.v_wt_r;  
    v_c_p( st_adptv.ca_v_k{lp_i+1}.v_k_idxs ) = v_c_p( st_adptv.ca_v_k{lp_i+1}.v_k_idxs ) + v_c_p_temp;  % Could probably be done cleaner (faster)
   
    
    % Now do eigenvector, if required
    if ( bp_calc_eigvec )
    
        % Weight the end of the main interval using weights from overlap LHS
        a_c_w_p( :, st_adptv.ca_overlap{lp_i}.v_k_idxs ) = a_c_w_p( :, st_adptv.ca_overlap{lp_i}.v_k_idxs ) .* st_adptv.ca_overlap{lp_i}.v_wt_l;

        % Get next interval
        a_c_w_p_temp = ca_dr{lp_i+1}.a_w_p;

        % Weight 
        a_c_w_p_temp( :, 1:Nwrt ) = a_c_w_p_temp( :, 1:Nwrt ) .* st_adptv.ca_overlap{lp_i}.v_wt_r;
        a_c_w_p( :, st_adptv.ca_v_k{lp_i+1}.v_k_idxs ) = a_c_w_p( :, st_adptv.ca_v_k{lp_i+1}.v_k_idxs ) + a_c_w_p_temp;
                
    end
    
end





end
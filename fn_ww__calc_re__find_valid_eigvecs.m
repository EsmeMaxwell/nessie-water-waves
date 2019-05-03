function [ c, v_w ] = fn_ww__calc_re__find_valid_eigvecs( v_zm, v_eig, a_w, U_min, U_max, ev_mode )
%fn_ww__calc_re__find_valid_eigvecs: Calc finds valid eigenvectors for Raleigh+free-surface problem
% 
%   [ c, v_w ] = fn_ww__calc_re__find_valid_eigvecs( v_zm, v_eig, a_w, U_min, U_max, ev_mode )
% 
% For the standard inviscid first-order Rayleigh equation with
% free-surface, the critical layer can cause problems. Here, we filter
% eigenvectors 
% 
% See also
%   fn_ww__calc_re__cl__red_c()



if ( nargin == 5 )
    ev_mode = 0;
end

if ( nargin < 5 || nargin > 6 )
    error( 'Wrong number of arguments' );
end
    

% TODO, in caller if we've hit here, don't need to redo the normalisation.

tol_imag = 1e-12;
tol_zero_bottom = 1e-12;
tol_surface_one = 1e-12;

% tol_zero_area_bottom = 0.1;
% zero_prop_bottom = 0.2;

Nev = size( a_w, 2 );
Nz = numel( v_zm );

interior_min_thrsh = -0.05;

area_factor = 1.05;

peak_min_idx = round( 0.05 * numel( v_zm ) );
peak_prominance = 0.04;




% Normalise first
a_w = fn_ww__util__normalise_eigenvectors( a_w, 0 );



% Eliminate those with nontrivial imag part (is rarely any use)
a_imag = abs( imag( a_w ) );
v_imag = sum( a_imag );
v_imag_idxs = find( v_imag > tol_imag );

v_eig( v_imag_idxs ) = [];
a_w( :, v_imag_idxs ) = [];


% Clean up a bit by ensuring everything is real
a_w = real( a_w );



% Eliminate those not near 0 at bottom
v_bottom = a_w( end, : );
v_bad_bottom_idxs = find( abs( v_bottom ) > tol_zero_bottom );

v_eig( v_bad_bottom_idxs ) = [];
a_w( :, v_bad_bottom_idxs ) = [];




% Remove lower branch (for now); or top branch depending on ev_mode
if ( 0 == ev_mode )
    v_lower_idxs = find( v_eig <= U_min );
    v_eig( v_lower_idxs ) = [];
    a_w( :, v_lower_idxs ) = [];
end



% Find those which aren't 1 at the surface
v_surf = a_w( 1, : );
v_surf_idxs = find( abs( v_surf - 1 ) > tol_surface_one );

v_eig( v_surf_idxs ) = [];
a_w( :, v_surf_idxs ) = [];






% Use MATLAB's peak finding code
if ( numel( v_eig ) > 1 )    
    v_crit_idxs = [];
    for lp_j=1:size(a_w,2)        
        [ v_peaks , v_peak_idxs ] = findpeaks( a_w(:,lp_j), 'MinPeakProminence', peak_prominance );
        if ( numel( v_peak_idxs ) > 0 && min( v_peak_idxs ) > peak_min_idx );
            v_crit_idxs(end+1) = lp_j;
        end            
    end
    
    % Only apply filter if we'll have at least one item left
    if ( numel( v_crit_idxs ) < numel( v_eig ) ) 
        v_eig( v_crit_idxs ) = [];
        a_w( :, v_crit_idxs ) = [];  
    end
end




% If we really need to, eliminate those which go sufificently -ve middle of
% domain, this is risky. Should also really chekc for prominance.
if ( numel( v_eig ) > 1 )
    v_min_val = min( a_w );
    v_neg_idxs = find( v_min_val < interior_min_thrsh );
    
    v_eig( v_neg_idxs ) = [];
    a_w( :, v_neg_idxs ) = [];    
end



% 
% % If we still have too many then we calcaulate areas and pick the least
% % area option
% if ( numel( v_eig ) > 1 )
%     v_area = trapz( abs( a_w ) );
%     [ min_area_val, min_area_idx ] = min( v_area ); 
% 
%     %v_bad_area_idxs = find( v_area > ( area_factor * min_area_val ) );
%     
%     %v_eig( v_bad_area_idxs ) = [];
%     %a_w( :, v_bad_area_idxs ) = [];
%     
%     v_eig = v_eig( min_area_idx );
%     a_w = a_w( :, min_area_idx );
% end


% TODO fix this.  Usually use the area thingy above but for the spectrum
% plots, we need to be able to choose lower or higher


% figure(1);
% size( v_zm )
% size( a_w )
% plot( v_zm, a_w );
% 


% Finally, if we still have too many, pick the dominant eigenvalue.
if ( numel( v_eig ) > 1 )
    if ( 0 == ev_mode ) 
        [ eig_max_val, eig_max_idx ] = max( v_eig );

        v_eig = v_eig( eig_max_idx );
        a_w = a_w( :, eig_max_idx );
    else
        [ eig_min_val, eig_min_idx ] = min( v_eig );

        v_eig = v_eig( eig_min_idx );
        a_w = a_w( :, eig_min_idx );
    end
end




if ( numel( v_eig ) < 1 )
    error( 'Woops, over-filtered the eigenvectors!' );
end

if ( numel( v_eig ) > 1 )
    v_eig
    plot( v_zm, a_w );
    error( 'Too many eigenvalues still left, cannot fix.  This really shouldn''t happen' );
end

c = v_eig(1);
v_w = a_w( :, 1 );


% % Eliminate those whose integral near bottom is over certain proporition of
% % its total integral
% cutoff_idx = floor( ( 1 - zero_prop_bottom ) * size( a_w, 1 ) )
% v_trapz_tot = trapz( abs( a_w ) )
% v_trapz_bottom = trapz( abs( a_w( cutoff_idx:end, : ) ) )
% v_bottom_area_prop = v_trapz_bottom ./ v_trapz_tot
% v_bad_bottom_area_idxs = find( v_bottom_area_prop  > tol_zero_area_bottom );





%v_invalid_idxs = union( v_imag_idxs , v_ev_negative_idxs );
%v_invalid_idxs = union( v_invalid_idxs, v_bad_bottom_idxs );



% v_invalid_idxs = union( v_invalid_idxs, v_crit_idxs );
% v_invalid_idxs = union( v_invalid_idxs, v_diff_negative_idxs )


% % Try to find eigenvectors which are from critical layers
% v_crit_idxs = [];
% for lp_j=1:size(a_w,2)
%     [ v_peaks , v_peak_idxs ] = findpeaks(a_w(:,lp_j),'MinPeakProminence',0.1);
%     if ( numel( v_peak_idxs ) > 0 && max ( v_peak_idxs ) < ( Nz - 1 ) )
%         v_crit_idxs(end+1) = lp_j;
%     end    
% end



% % Eliminate those which are not monotonically increasing
% a_diff = diff( -a_w, 1, 1 );
% v_diff_min = min( a_diff );
% v_diff_negative_idxs = find( v_diff_min < 0 )


% hold on;
% for lp_j=1:size(a_w_valid,2)
%     [ v_peaks , v_peak_idxs ] = findpeaks(a_w_valid(:,lp_j),'MinPeakProminence',0.1);
%     scatter( v_zm( v_peak_idxs ), v_peaks );
% end
% hold off;



        
%         if ( lp_k == 37 )
%     %         a_w = fn_ww__util__normalise_eigenvectors( a_w, 1 );
%     %         a_w = abs( a_w );
%     %         plot( st_Dn.v_zm, a_w )
%             [ v_r, v_cond ] = fn_ww__calc_re__resolvent_c( a_A2, a_A1, a_A0, v_eig, 1, 1, 1 );
%             semilogy( real( v_eig ), v_cond, 'r*' );
%             v_c_sp = linspace( -1.2, 0.2, 10000 );
%             [ v_r_sp, v_cond_sp ] = fn_ww__calc_re__resolvent_c( a_A2, a_A1, a_A0, v_c_sp, 1, 1, 1 );
%             hold on;
%             semilogy( v_c_sp, v_cond_sp, 'b' );
%             hold off;       
%             error();        
%         end


end
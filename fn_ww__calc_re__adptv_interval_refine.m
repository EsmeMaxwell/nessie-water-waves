function [ st_adptv ] = fn_ww__calc_re__adptv_interval_refine( v_k, st_p )
%fn_ww__calc_re__cl_adptv__red_c: Calc adaptive intervals for dispersion relation calcs
% 
% 
% 
% 
% TODO need to be careful with normalisation, cf. pf

tau = st_p.fp_adtpv_depth_tau;

k_min = min( v_k );
k_max = max( v_k );

ltol = -log( st_p.fp_adtpv_depth_tol );
hmin = st_p.fp_adtpv_depth_hmin;
hmax = st_p.fp_adtpv_depth_hmax;

% If we only need one interval, do that; if not, then we need to setup the
% first interval in a nondefault manner so cannot loop over it.
if ( ltol > k_max )
    v_intervals = [ 0 ltol / hmin ];
    v_optimal_h = [ st_p.h ];
else
    v_intervals = [ 0 ltol/ hmax; ...
                    ltol / 1.0, ltol / ( hmin * hmax ) ];
    v_optimal_h = [ st_p.h; min( [ st_p.h, 0.5 * ( 1 + hmin * hmax ) ] ) ];
end


% Generate overlapping intervals and calculate optimal h
lp_ctr = 1;
while ( v_intervals( end, 2 ) < k_max )
    lp_ctr = lp_ctr + 1;    
    v_intervals = [ v_intervals; ltol / hmin^(lp_ctr-1), ltol/ ( hmin^(lp_ctr) * hmax )  ];
    v_optimal_h = [ v_optimal_h; 0.5 * ( hmin^(lp_ctr-1) + ( hmin^(lp_ctr) * hmax ) ) ];    % TODO for safety, should this be min{h, ...}?
end

% Clean-up
if ( v_intervals( end, 2 ) > k_max )
    v_intervals( end, 2 ) = k_max;
end
Nintervals = size( v_intervals, 1 );



% Calculate overlap regions
Noverlaps = max( [ 0, Nintervals - 1 ] );
v_overlap = zeros( Noverlaps, 2 );
for lp_i=2:Nintervals
    v_overlap(lp_i-1,:) = [ v_intervals(lp_i,1), v_intervals(lp_i-1,2 ) ];
end



% Now actually assign k points into intervals.
ca_v_k = cell( 1, Nintervals );
for lp_i=1:Nintervals
    v_k_idxs = find( v_k >= v_intervals(lp_i,1) & v_k <= v_intervals(lp_i,2) );
    ca_v_k{lp_i}.v_k = v_k( v_k_idxs );
    ca_v_k{lp_i}.v_k_idxs = v_k_idxs;    
end

% Now assign into overlap intervals
if ( Noverlaps > 0 )
    
    %Noverlaps = Nintervals - 1;
    ca_overlap = cell( 1, Noverlaps );
    
    for lp_i=1:Noverlaps
        v_k_idxs = find( v_k >= v_overlap(lp_i,1) & v_k <= v_overlap(lp_i,2) );
        ca_overlap{lp_i}.v_k = v_k( v_k_idxs );
        ca_overlap{lp_i}.v_k_idxs = v_k_idxs;
       
        %numel( v_k_idxs )
        assert( numel( v_k_idxs ) >= st_p.ip_adtpv_min_numel_overlap, 'Insufficient elements in overlap region' );            
    end
   
end

    

% Now calculate parition of unity weights and apply for each of our overlap
% intervals

% First define the phi function (see Aiton & Discoll 2017). Must ensure
% compact in [-1,1]. Also declare the weight functions, etc.
fn_ind = @(x) ( abs(x) < 1 );
fn_phi = @(v_x) fn_ind(v_x) .* exp( 1 - 1 ./ ( 1 - v_x.^2 ) );
fn_phi_l = @(v_x,tau) fn_phi( ( v_x + 1 ) ./ ( 1 + tau ) );
fn_phi_r = @(v_x,tau) fn_phi( ( v_x - 1 ) ./ ( 1 + tau ) );
fn_wt_l = @(v_x,tau) fn_phi_l(v_x,tau) ./ ( fn_phi_l(v_x,tau) + fn_phi_r(v_x,tau) );
fn_wt_r = @(v_x,tau) fn_phi_r(v_x,tau) ./ ( fn_phi_l(v_x,tau) + fn_phi_r(v_x,tau) );


for lp_i=1:Noverlaps
   
    % Map to [-1,1].
    v_overlap_z0 = fn_ww__setup__lin_conv_zm_to_z0( ca_overlap{lp_i}.v_k, v_overlap(lp_i,1), v_overlap(lp_i,2)  );
    
    % Calculate weights
    ca_overlap{lp_i}.v_wt_l = fn_wt_l( v_overlap_z0, tau );
    ca_overlap{lp_i}.v_wt_r = fn_wt_r( v_overlap_z0, tau );
    
end


st_adptv = struct;
st_adptv.v_k = v_k;
st_adptv.v_intervals = v_intervals;
st_adptv.ca_v_k = ca_v_k;
st_adptv.v_optimal_h = v_optimal_h;
st_adptv.v_overlap = v_overlap;
st_adptv.ca_overlap = ca_overlap;
st_adptv.Nintervals = Nintervals;
st_adptv.Noverlaps = Noverlaps;


% plot( ca_overlap{1}.v_k, ca_overlap{1}.v_wt_l, 'm--', ca_overlap{1}.v_k, ca_overlap{1}.v_wt_r, 'g-.', ca_overlap{1}.v_k, ca_overlap{1}.v_wt_l + ca_overlap{1}.v_wt_r, 'b' )

end
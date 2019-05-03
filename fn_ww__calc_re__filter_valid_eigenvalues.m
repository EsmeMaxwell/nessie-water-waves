function [ v_eig, a_w, v_valid_idxs ] = fn_ww__calc_re__filter_valid_eigenvalues( v_eig, a_w, st_p )
%fn_ww__calc_re__filter_valid_eigenvalues: Calc filters suitable c eigenvalues
%
%   [ v_eig, a_w, v_valid_idxs ] = fn_ww__calc_re__filter_valid_eigenvalues( v_eig, a_w, st_p )
% 
% Filters sensible eigenvalues when c is eigenvalue.
%
% See also
%   fn_ww__calc_re__cl__red_c()


v_nonfinite_idxs = find( ~isfinite( v_eig ) ).';
v_unsafe_idxs = find( abs( v_eig ) > st_p.fp_safety_thrsh ).';
v_imag_idxs = find( abs( imag( v_eig ) ) > st_p.fp_imag_thrsh ).';

v_bad_idxs = union( v_nonfinite_idxs, v_unsafe_idxs );
v_bad_idxs = union( v_bad_idxs, v_imag_idxs );

v_valid_idxs = 1:numel(v_eig);
v_valid_idxs = setdiff( v_valid_idxs, v_bad_idxs );


v_eig = v_eig( v_valid_idxs );
a_w = a_w( :, v_valid_idxs );

    

end
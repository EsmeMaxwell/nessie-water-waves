function [ st_r_shear ] = fn_ww__setup__create_shear_r_st__fn( st_fn_shear, st_p )
%fn_ww__setup__create_shear_r: Setup return a reduced shear rofile struct from st_fn_shear
%
%   [ st_r_shear ] = fn_ww__setup__create_shear_r( st_fn_shear )
%
% Almost trivial, just stores the shear fns in another struct along with
% Umin, Umax results.
%
% INPUT
%   st_fn_shear : Struct of shear profile functions
%
% OUTPUT
%   st_r_shear : Struct containing shear functions and Umin,Umax inf
%
% TAGS: CORE, SISCPFLIB, WWERRINSHEAR
%
% See also
%   fn_ww__setup__create_shear_r_st__vec(),
%   fn_ww__util__get_r_shear_data(),
%   fn_ww__setup__param_std__re_cl()


% Copy shear functions
st_r_shear = struct;
st_r_shear.st_fn_shear = st_fn_shear;

if ( st_p.bp_Uminmax_precalc )    
    % Calculate function U_max, U_min
    [ st_crit ] = fn_ww__util__find_U_min_max__fn( st_fn_shear, st_p );
    st_r_shear.st_crit = st_crit;
end


end
function [ st_vec ] = fn_ww__setup__create_st_vec( v_zm, v_zp )
%fn_ww__setup__create_st_vec: Setup create st_vec from one or two vectors, almost redundant
%
%   [ st_vec ] = fn_ww__setup__create_st_vec( v_zm, v_zp )
%

st_vec = struct;
st_vec.v_zm = v_zm;

if ( nargin > 1 )
    st_vec.v_zp = v_zp;
end


end
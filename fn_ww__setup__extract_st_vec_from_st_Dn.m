function [ st_vec ] = fn_ww__setup__extract_st_vec_from_st_Dn( st_Dn )
%fn_ww__setup__extract_st_vec_from_st_Dn: Setup create st_vec from st_Dn
%
%   [ st_vec ] = fn_ww__setup__extract_st_vec_from_st_Dn( st_Dn )
%

st_vec = struct;
st_vec.v_zm = st_Dn.v_zm;

if ( isfield( st_Dn, 'v_zp' ) ) % In fairness, this should _always_ be true
    st_vec.v_zp = st_Dn.v_zp;
end


end
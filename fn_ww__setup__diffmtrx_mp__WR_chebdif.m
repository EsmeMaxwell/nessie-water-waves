function [ st_Dn ] = fn_ww__setup__diffmtrx_mp__WR_chebdif( N, mp_digits )
%fn_ww__setup__diffmtrx_mp__WR_chebdif: Generate diff matrices using WR chebdif in mp
%
%   [ st_Dn ] = fn_ww__setup__diffmtrx_mp__WR_chebdif( N, b_nst )
% 
% High precision version to be used with Advanpix library. Does the same as
% FN_WW__SETUP__DIFFMTRX__WR_CHEBDIF() but using appropriate mp() syntax.
%
% Constructs the first and second order differentiation matrix, of size
% (N+1), using chebfif from the Weideman-Reddy DM suite. The negative sum
% trick can be applied by setting the b_nst flag nonzero. Also returns the
% corresponding z vector of second-kind Chebyshev collocation points on
% [-1,1].
% 
% INPUT
%
% N : Size, actually will return a matrix and vector of size (N+1)
% 
% b_nst : Flag for whether to apply negative sum trick. Anything nonzero
% will enable.
% 
% OUTPUT
%
% st_Dn : struct containing matrices a_D, a_D2, and v_z0.
%
% See also
%   fn_ww__setup__diffmtrx__WR_chebdif(),
%   fn_ww__setup__diffmtrx__WR_poldif(),
%

if ( nargin > 1 )
    mp.Digits( mp_digits );
else
    mp_digits = 34;
    mp.Digits( mp_digits );
end
    
[ v_z0, a_DM] = fn_ww__ext__diffmtrx_mp__WR_chebdif(N+1, 2);
a_D = a_DM(:,:,1);
a_D2 = a_DM(:,:,2);

st_Dn = struct;
st_Dn.v_z0 = v_z0;
st_Dn.a_D = a_D;
st_Dn.a_D2 = a_D2;
st_Dn.mp_digits = mp_digits;
st_Dn.N = N;

end
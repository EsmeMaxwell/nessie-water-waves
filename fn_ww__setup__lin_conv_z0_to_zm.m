function [ v_zm ] = fn_ww__setup__lin_conv_z0_to_zm( v_z0, a, b )
%fn_ww__setup__lin_conv_z0_to_zm: Setup map collocation to nondimensional coordinates
%
%   [ v_zm ] = fn_ww__setup__lin_conv_z0_to_zm( v_z0, a, b )
%
% Simply maps the supplied collocation coordinates to the specified
% physical interval, [a,b].
% 
% INPUT
%   v_z0 : Vector of collocation points
%   a, b : Physical interval to map to
% 
% OUTPUT
%   v_zm : Vector of the coordinates in the physical interval
% 
% 
% See also
%   fn_ww__setup__lin_conv_zm_to_z0()


v_zm =  0.5 * ( v_z0 * ( b - a ) + ( b + a ) );

end
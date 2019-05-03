function [ v_z0 ] = fn_ww__setup__lin_conv_zm_to_z0( v_zm, a, b )
%fn_ww__setup__lin_conv_zm_to_z0: Setup map nondimeinsional coordinates back to collocation domain
%
%   [ v_z0 ] = fn_ww__setup__lin_conv_zm_to_z0( v_zm, a, b )
%
% Sipmly maps the supplied nondimeinsional coordinates to the collocation
% interval, [-1,1]
% 
% INPUT
%   v_zm : Vector of the coordinates in the physical interval
%   a, b : Physical interval used
% 
% OUTPUT
%   v_z0 : Vector of collocation points in the [-1,1] interval
% 
% See also
%   fn_ww__setup__lin_conv_z0_to_zm()

v_z0 = ( 2 * v_zm - ( b + a ) ) / ( b - a );

end
function [ ver_maj, ver_min, ver_serial, s_notice ] = fn_ww__setup__version()
%fn_ww__setup__version: Setup return current version
% 
%   [ ver_maj, ver_min, ver_serial, s_notice ] = fn_ww__setup__version()
% 
% Returns the current version information.
%
% OUTPUT
%   ver_maj : Major version number as integer
%   ver_min : Minor version number as integer
%   ver_serial : Increment within release
%   s_notice : Type of revision
% 

    ver_maj = 0;
    ver_min = 1;
    ver_serial = 15;
    s_notice = 'dev-restricted';

    fprintf( 'Version %d.%d.%d.%s\n', ver_maj, ver_min, ver_serial, s_notice );
    
end
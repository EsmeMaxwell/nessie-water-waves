function [ st_p ] = fn_ww__setup__merge_parameters( st_p, st_p_ctl )
%fn_ww__setup__merge_parameters: Setup merge new parameters into existing parameter set
%
%   [ st_p ] = fn_ww__setup__merge_parameters( st_p, st_p_ctl )
%
% Merge parameters from st_p_ctl into st_p.
%
%
% INPUT
%   st_p : Parameter set to merge into
%   st_p_ctl : Parameter set with new parmeters
% 
% OUTPUT
%   st_p : Parameter set after merge
%
%
% See also
%   fn_ww__setup__param_std__re_cl(),
%   fn_ww__setup__param_ctl__re_cl(),
%   fn_ww__setup__param_ctl__re_cl__exp()

s_names = fieldnames( st_p_ctl );

for lp_i=1:numel( s_names )

    % Recursion check
    if ( isfield( st_p, s_names{lp_i} ) && isstruct( st_p.(s_names{lp_i}) ) )
        
        % We can only submerge in other structs
        if ( ~isstruct( st_p_ctl.(s_names{lp_i}) ) )
            error( 'Must overwrite struct with another struct.');
        end
        
        [ st_sub ] = fn_ww__setup__merge_parameters( st_p.(s_names{lp_i}), st_p_ctl.(s_names{lp_i}) );
        st_p.(s_names{lp_i}) = st_sub;
    
    else
    
        st_p.(s_names{lp_i}) = st_p_ctl.(s_names{lp_i});
        
    end
    
end


end
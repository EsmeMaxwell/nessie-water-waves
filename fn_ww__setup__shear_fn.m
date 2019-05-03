function [ st_fn_shear ] = fn_ww__setup__shear_fn( profile_id, ca_args )
%FN_WW__SETUP__SHEAR_FN: Choose from a set of standard shear profiles
%
% [ st_fn_shear ] = fn_ww__setup__shear_fn( profile_id, ca_args )
%
% Returns a struct containing functions fn_U(z), fn_dU(z), fn_ddU(z) along
% with parameters which generate specific shear functions. A selection of
% standard profiles are supplied and the user can alter the parameters
% using ca_args.
%
% All the shear profiles here are on a physical domain and assume the z
% interval is [-1,0]. If you want to use a different interval, you'll need
% to adjust manually.
%
% INPUT
%
% profile_id : An integer to select from the available shear profiles. The
% table below enumerates,
%   10=Translated cos(z^2) profile,
%       'csp_b * ( cos( csp_a * (z-csp_z0).^2 + csp_z1 ) ) + csp_c'   
%   11=Translated -cos(z^2) profile
%   20=Translated exponential profile
%       'esp_b * ( exp( esp_lambda * z ) ) + esp_c0'
%
%
% ca_args : Argument list, must be in form of { 'itemA', valueA, 'itemB',
% valueB }. Use variable names listed above.
%
% OUTPUT
%
% st_fn_shear : struct containing the functions and parameters.
% 
% See also 
%   FN_WW__SETUP__SHEAR_FN_TO_VEC(),
%   FN_WW__SETUP__LIN_CONV_Z0_TO_ZP(), 
%   FN_WW__SETUP__LIN_CONV_ZP_TO_Z0()
%

% TODO self-consistency check for derivatives here.


FrS = 0.1;
U0 = 0.0;

csp_a = 4*pi/7;
csp_b = 0.2;
csp_c = 0.2;
csp_t = 0;
csp_z0 = 0.3;
csp_z1 = 0;

esp_lambda = 1;
esp_b = 1;
esp_c0 = 0;


st_fn_shear = struct;
st_fn_shear.b_lin_shear = 0;

% Process input arguments
for lp_k=1:2:length(ca_args)
    switch(ca_args{lp_k})
        case 'FrS'
            FrS = ca_args{lp_k+1};
        case 'U0'
            U0 = ca_args{lp_k+1};
        case 'csp_a'
            csp_a = ca_args{lp_k+1};
        case 'csp_b'
            csp_b = ca_args{lp_k+1};
        case 'csp_c'
            csp_c = ca_args{lp_k+1};            
        case 'csp_t'
            csp_t = ca_args{lp_k+1};
        case 'csp_z0'
            csp_z0 = ca_args{lp_k+1};
        case 'csp_z1'
            csp_z1 = ca_args{lp_k+1};
        case 'esp_lambda'            
            esp_lambda = ca_args{lp_k+1};            
        case 'esp_b'
            esp_b = ca_args{lp_k+1};            
        case 'esp_c0'
            esp_c0 = ca_args{lp_k+1};                        
        otherwise
            error( 'Invalid parameter selection setting up shear profile.' );
    end
end



switch profile_id        
    case 10
        % Translated cos(z^2) profile (this IS different to case 3, fixed an error)
        st_fn_shear.fn_U = @(z) csp_b * ( cos( csp_a * (z-csp_z0).^2 + csp_z1 ) ) + csp_c;
        st_fn_shear.fn_dU = @(z) csp_b * ( -2 * csp_a * (z-csp_z0) .* sin( csp_a * (z-csp_z0).^2 + csp_z1 ) );
        st_fn_shear.fn_ddU = @(z) csp_b * ( -4 * csp_a^2 * (z-csp_z0).^2 .* cos( csp_a * (z-csp_z0).^2 + csp_z1 ) -2 * csp_a * sin( csp_a * (z-csp_z0).^2 + csp_z1 ) );
        st_fn_shear.csp_a = csp_a;
        st_fn_shear.csp_b = csp_b;
        st_fn_shear.csp_c = csp_c;        
        st_fn_shear.csp_z0 = csp_z0;
        st_fn_shear.csp_z1 = csp_z1;

        
    case 11
        % Translated -cos(z^2) profile (this IS different to case 3, fixed an error)
        st_fn_shear.fn_U = @(z) -csp_b * ( cos( csp_a * (z-csp_z0).^2 + csp_z1 ) ) + csp_c;
        st_fn_shear.fn_dU = @(z) -csp_b * ( -2 * csp_a * (z-csp_z0) .* sin( csp_a * (z-csp_z0).^2 + csp_z1 ) );
        st_fn_shear.fn_ddU = @(z) -csp_b * ( -4 * csp_a^2 * (z-csp_z0).^2 .* cos( csp_a * (z-csp_z0).^2 + csp_z1 ) -2 * csp_a * sin( csp_a * (z-csp_z0).^2 + csp_z1 ) );
        st_fn_shear.csp_a = csp_a;
        st_fn_shear.csp_b = csp_b;
        st_fn_shear.csp_c = csp_c;        
        st_fn_shear.csp_z0 = csp_z0;
        st_fn_shear.csp_z1 = csp_z1;        
        
    case 20
        % Translated exponential profile
        st_fn_shear.fn_U = @(z) esp_b * ( exp( esp_lambda * z ) ) + esp_c0;
        st_fn_shear.fn_dU = @(z) esp_b * ( esp_lambda * exp( esp_lambda * z ) );
        st_fn_shear.fn_ddU = @(z) esp_b * ( esp_lambda^2 * exp( esp_lambda * z ) );
        st_fn_shear.esp_lambda = esp_lambda;
        st_fn_shear.esp_b = esp_b;        
        st_fn_shear.esp_c0 = esp_c0;
        
        
        
%     case 20
%         % Log profile
%         st_CL_fn_shear.fn_U = @(z) 
%         st_CL_fn_shear.fn_dU = @(z) 
%         st_CL_fn_shear.fn_ddU = @(z)
%         
        
    otherwise
        error( 'No valid shear profile selected' );
end





end
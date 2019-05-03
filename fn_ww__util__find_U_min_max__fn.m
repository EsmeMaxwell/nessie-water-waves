function [ st_crit ] = fn_ww__util__find_U_min_max__fn( st_fn_shear, st_p )
%fn_ww__util__find_U_min_max: Util find Umin and Umax of fn shear profile
% 
%   [ st_crit ] = fn_ww__util__find_U_min_max__fn( st_fn_shear, st_p )
% 
% 
% 


st_crit = struct;


fn_U_min = @(z) (st_fn_shear.fn_U(z));
fn_U_max = @(z) -(st_fn_shear.fn_U(z));

x0 = ( st_p.a + st_p.b ) / 2;
ob_opts = optimoptions( @fmincon, 'Algorithm', 'sqp', 'Display', 'off', ... 
        'OptimalityTolerance', 1e-11, ...
        'StepTolerance', 1e-11, ...
        'ConstraintTolerance', 1e-11, ...
        'MaxFunctionEvaluations', 300, ...
        'MaxIterations', 300 );
ob_gs = GlobalSearch( 'Display', 'off', ...
                        'XTolerance', 1e-10, ...
                        'FunctionTolerance', 1e-10, ...
                        'NumStageOnePoints', 15, ...
                        'NumTrialPoints', 50 );

% Find min
st_problem_min = createOptimProblem( 'fmincon', 'objective', fn_U_min, 'x0', x0, 'lb', st_p.a, 'ub', st_p.b, 'options', ob_opts );
[ z_min, U_min ] = run( ob_gs, st_problem_min );

% Find max
st_problem_max = createOptimProblem( 'fmincon', 'objective', fn_U_max, 'x0', x0, 'lb', st_p.a, 'ub', st_p.b, 'options', ob_opts );
[ z_max, U_max ] = run( ob_gs, st_problem_max );
U_max = -1 * U_max;

% Collate
st_crit.z_min = z_min;
st_crit.z_max = z_max;
st_crit.U_min = U_min;
st_crit.U_max = U_max;

% If reuqired, calculate physicsl coords
if ( isfield( st_fn_shear, 'fn_phy_U' ) )
    st_crit.z_phy_min = st_crit.z_min * st_p.phy_h;
    st_crit.z_phy_max = st_crit.z_max * st_p.phy_h;
    st_crit.U_phy_min = st_fn_shear.fn_phy_U( st_crit.z_phy_min );
    st_crit.U_phy_max = st_fn_shear.fn_phy_U( st_crit.z_phy_max );
end




% 
% 
% % Try to determine whether we've got a numeric or function for shear
% % profile
% if ( isnumeric( var_U ) && numel( var_U ) > 1 ) 
%     
%     v_U = var_U;        
%     [ U_min ] = min( v_U );
%     [ U_max ] = max( v_U );
%     
%     % Collate and calculate physical
%     st_crit.U_min = U_min;
%     st_crit.U_max = U_min;
%     st_crit.U_phy_min = st_fn_shear.fn_phy_U( st_crit.z_phy_min );
%     st_crit.U_phy_max = st_fn_shear.fn_phy_U( st_crit.z_phy_max );
%     
% elseif ( isa( var_U, 'function_handle' ) )   
%     
%     fn_U_min = @(z) var_U(z);
%     fn_U_max = @(z) -var_U(z);
%     
%     x0 = ( st_p.a + st_p.b ) / 2;
%     ob_opts = optimoptions( @fmincon, 'Algorithm', 'sqp', 'Display', 'off', ... 
%             'OptimalityTolerance', 1e-12, ...
%             'StepTolerance', 1e-12, ...
%             'ConstraintTolerance', 1e-12, ...
%             'MaxFunctionEvaluations', 500, ...
%             'MaxIterations', 500 );
%     ob_gs = GlobalSearch( 'Display', 'off', 'XTolerance', 1e-10, 'FunctionTolerance', 1e-10 );    
%     
%     % Find min
%     st_problem_min = createOptimProblem( 'fmincon', 'objective', fn_U_min, 'x0', x0, 'lb', st_p.a, 'ub', st_p.b, 'options', ob_opts );
%     [ z_min, U_min ] = run( ob_gs, st_problem_min );
%     
%     % Find max
%     st_problem_max = createOptimProblem( 'fmincon', 'objective', fn_U_max, 'x0', x0, 'lb', st_p.a, 'ub', st_p.b, 'options', ob_opts );
%     [ z_max, U_max ] = run( ob_gs, st_problem_max );
%     U_max = -1 * U_max;
% 
%     % Collate and calculate physical
%     st_crit.z_min = z_min;
%     st_crit.z_max = z_max;
%     st_crit.U_min = U_min;
%     st_crit.U_max = U_min;
%     st_crit.z_phy_min = st_crit.z_min * st_p.phy_h;
%     st_crit.z_phy_max = st_crit.z_max * st_p.phy_h;
%     st_crit.U_phy_min = st_fn_shear.fn_phy_U( st_crit.z_phy_min );
%     st_crit.U_phy_max = st_fn_shear.fn_phy_U( st_crit.z_phy_max );
%             
% else
%     error( 'Must supply either a vector or function handle.' );
% end
% 


end
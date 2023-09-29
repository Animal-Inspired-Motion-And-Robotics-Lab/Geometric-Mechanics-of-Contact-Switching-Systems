function Mp = Legged_dissipation_metric(geometry,physics,shapeparams)
% Calculate the dissipation power metric for a set of curvature bases

% Identify what kind of system is being calculated, and use this to specify 
% how the dissipation metric should be generated
switch geometry.type
    
    % Check if this is a single disk-legged type
    case 'n-disk-legged'

        switch physics.type
    
            case 'isotropic_friction'
                physics_function = @viscTranslationalFriction_dissipation_metric;
        
            case 'anisotropic_friction'
                error('ERROR: Viscous friction anisotropy is not supported at the moment!');

            case 'inertial'
                error('ERROR: Inertial metric is not supported at the moment!');

            case 'rate_limited'
                physics_function = @Rate_Limitation_metric;
                
        end
    
    % If it is a multdisk legged system --  legs don't have the same
    % center of rotation (the system can be discretely compliant in the 
    % CLARI-sense or still have a rigid body with offset links
    case 'n-multidisk-legged'

        error('ERROR: Multidisk systems are not supported at the moment!');

%         switch physics.type
%     
%             case 'isotropic'
%                 physics_function = @MultiDisk_viscTransFric_metric;
%                 error('ERROR: Multidisk systems are not supported at the moment!');
%         
%             case 'anisotropic'
%                 error('ERROR: Multidisk systems are not supported at the moment!');
% 
%             case 'inertial'
%                 error('ERROR: Multidisk systems are not supported at the moment!');
%                 
%         end

end

% Call the physics function identified for the system
Mp = physics_function(geometry,physics,shapeparams);

end


% 	% Specified integration limits
% 	int_limit = geometry.length*[-0.5 0.5];
% 	
% 	% Define the tangential, lateral drag matrix for unit/double drag
% 	drag_matrix = [1 0; 0 physics.drag_ratio]*physics.drag_coefficient;
% 
% 	% Get the backbone locus, Jacobian, and Local Connection functions
% 	[A, h, J] = LowRE_local_connection(geometry,physics,shapeparams);
% 
% 	% Integrate along the body for the power metric
% 	Mp_sol = ode45(@(s,Mp) dMetric(s,Mp,A,h(s),J(s),drag_matrix),int_limit,zeros(length(shapeparams)));%zeros(length(shapeparams)^2,1));
% 
%     % Form the components of the metric into the standard NxN symmetric matrix form
% 	Mp = reshape(deval(Mp_sol,int_limit(end)),length(shapeparams),[]);
% 
% end
% 
% function dMp = dMetric(s,Mp,A,h,J,drag_matrix) %#ok<INUSL>
% % Calculate the local contribution to the power metric
% 
%     % Calculate the jacobian from shape variables to body velocity of this
%     % point on the system:
% 	
%     % Body velocity for a point on the robot is sum of right-translated
%     % system body velocity and velocity within body frame, left translated
%     % into local coordinates
%     A_point = -(TgLginv(h) * (-TeRg(h)*A + J));
%     
%     % the xy jacobian from shape variables to body velocity of points on the system is the
%     % first two rows of the local connection.
%     localJ = -A_point(1:2,:);
% 
%     % Contribution to the metric from this point is its local drag matrix,
%     % pulled back into joint coordinates by its Jacobian.
% 	dMp = localJ.'*drag_matrix*localJ;
% 	
%     % Convert metric contribution into column vector for ODE integration
% 	dMp = dMp(:);
% 	
% end
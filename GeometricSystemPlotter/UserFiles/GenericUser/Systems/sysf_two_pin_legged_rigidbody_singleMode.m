function output = sysf_two_pin_legged_rigidbody_singleMode(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Pin-Legged Rigid System: 2-legged -- V1.0 Bad'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/RigidDiskLegs/',...
                'Physics/LeggedSystems/'});
            
		case 'initialize'

            %%%%%%
            % Ankle limit amplitude for the legged robot:
            ank = pi/2;

            % Leg-to-link length fraction:
            a = 1;

            % Set the body-length of the robot:
            bl = 1;

            % Define system geometry
            s.geometry.type = 'n-disk-legged'; % Other types NOT SUPPORTED RN -- 'n-multidisk-legged'
            s.geometry.linklengths = bl*[a/2 1 a/2]; % body - middle link
            s.geometry.baseframe = 'center';
            s.geometry.length =...
                sum(s.geometry.linklengths)/(numel(s.geometry.linklengths)-1); 
            % let's get the body-link to have length 1
            s.geometry.link_shape = {'ellipse','ellipse','ellipse'};
                st_legs = struct('aspect_ratio',0.1);
                st_body = struct('aspect_ratio',0.5);
            s.geometry.link_shape_parameters = {st_legs,st_body,st_legs};
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = ank*[-1  0  1];

            %%%
            %%%%%%
            % Here, we directly move on to describing the system physics,
            % but in a more general case, we would like to evaluate the
            % actuator transmission model to obtain friction guarantees and
            % then we can set the 'static-friction-coefficients' based on
            % use case, etc. We need an actuator2physics function where we 
            % describe the system.

            % Define system physics
            s.physics.st_friction_coeff = [0, nan, inf]; % 0 --for slipping 
                                                         % legs
                                                         % nan -- for
                                                         % rigid-body non
                                                         % contacting links
                                                         % inf -- for
                                                         % pinned legs
            % Other cases for the system - [inf, nan, 0], [0, nan, 0]
            % The second case would be interesting since that would purely
            % be frictional motion.

            s.physics.coulomb_friction_coeff = [];
            s.physics.viscous_friction_coeff = 1;
            s.physics.type = 'isotropic_friction';
            % Other types include (NOT SUPPORT RN) ~~~~~~~~~~~~~~~~~~~~~~~~
            % 'anisotropic_friction' 'inertial'
           
            % Functional Local connection
                    % in physics
            s.A = @(alpha1,alpha2) Legged_local_connection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2]);                        % Joint angles
                    % in physics
            s.metric = @(alpha1,alpha2) Legged_dissipation_metric(...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2]);                        % Joint angles

            % No dissipation metric for now. Need to incorporate this when
            % we start slipping.
     
			%%%
            
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1]*ank; % let's say the leg-amp limit is 45degs

			%densities for various operations
			s.density.vector = [21 21 ]; %density to display vector field
			s.density.scalar = [51 51 ]; %density to display scalar functions
			s.density.eval = [31 31 ];   %density for function evaluations
            s.density.metric_eval = [11 11]; %density for metric evaluation
            s.density.finite_element=31;
            s.density.coriolis_eval = [31 31];
            s.density.mass_eval = [31 31]; % density for mass matrix evaluation

			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*ank;
            s.tic_locs.y = [-1 0 1]*ank;

            % Set system type variable for gait optimization
            s.system_type = 'drag';
            
			%%%%
			%Save the system properties
			output = s;

	end

end


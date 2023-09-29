function output = sysf_two_pin_legged_rigidbody_SmoothModeChange(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Pin-Legged Rigid System: 2-legged, Smooth Contact Switch'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Geometry/RigidDiskLegs/',...
                'Physics/LeggedSystems/'});
            
		case 'initialize'

            %%%%%%
            % Ankle limit amplitude for the legged robot:
            ank = round(pi/2,1);

            % Contact map domain limits:
            % OLD DOMAIN: [0, 1]
            con_lim = [-1, 1]; % This won't change based on the function.
            dom_margin = 10;  % in percentage give how much larger the input
                             % space needs to be incomparison to the
                             % contact domain.
            dom_thresh = [con_lim(1)-dom_margin*diff(con_lim)/100,...
                con_lim(2)+dom_margin*diff(con_lim)/100];
            
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
            s.physics.st_friction_coeff = [0, nan, inf;
                                           inf, nan, 0]; % 0 --for slipping 
                                                         % legs
                                                         % nan -- for
                                                         % rigid-body non
                                                         % contacting links
                                                         % inf -- for
                                                         % pinned legs
            % We have specified the two modes we are interested in since we
            % shall try to move smoothly from one mode to the other as the
            % legs swap their adhesion strengths.

            s.physics.coulomb_friction_coeff = [];
            s.physics.viscous_friction_coeff = [];
            s.physics.type = 'rate_limited';
            
            % Let's introduce another array that holds the metric main
            % diagonal entries in increasing order
%             s.physics.mDiag = [0.1, 1];
%             % If you want them to weighted equally:
            s.physics.mDiag = [];

            % Other types include (NOT SUPPORT RN) ~~~~~~~~~~~~~~~~~~~~~~~~
            % 'anisotropic_friction' 'inertial'
           
            % Functional Local connection
                    % in physics
            s.A = @(alpha,c) Legged_local_connection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha,c]);                        % Joint angles
                    % in physics
            s.metric = @(alpha,c) Legged_dissipation_metric(...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha,c]);                        % Joint angles and contact

            % No dissipation metric for now. Need to incorporate this when
            % we start slipping.
     
			%%%
            
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-ank,ank,dom_thresh(1),dom_thresh(2)];
            
			%densities for various operations
			s.density.vector = [21 21]; %density to display vector field
			s.density.scalar = [51 51 ]; %density to display scalar functions
			s.density.eval = [31 31 ];   %density for function evaluations
            s.density.metric_eval = [11 11]; %density for metric evaluation
            s.density.finite_element=31;
            s.density.coriolis_eval = [31 31];
            s.density.mass_eval = [31 31]; % density for mass matrix evaluation

			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*ank;
            s.tic_locs.y = [con_lim(1) 0 con_lim(2)];
            
            % Set system type variable for gait optimization
            s.system_type = 'drag'; % need to think about the choices
            
			%%%%
			%Save the system properties
			output = s;

	end

end
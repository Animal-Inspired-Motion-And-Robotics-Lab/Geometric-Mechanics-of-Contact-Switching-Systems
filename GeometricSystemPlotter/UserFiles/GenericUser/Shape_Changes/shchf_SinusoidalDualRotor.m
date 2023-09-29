function output = shchf_SinusoidalDualRotor(input_mode,pathnames)

	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end	
	
	switch input_mode
		
		case 'name'
			
			output = 'Double Rotor DS Multi-disk: Circle, Amp = 0.25, Offset = [0.3,-0.3]';
			
		case 'dependency'
			
			output.dependency = {};
			
		case 'initialize'

			%%%%
			% Filename to save to

			%%
			%Path definitions

			%path definition
			p.phi_def = @strokedef;
			
			%marker locations
			p.phi_marker = [];
			
			%arrows to plot
			p.phi_arrows = 2;

			%time to run path
			p.time_def = [0 2*pi];

			p.cBVI_method = 'simple';

			%path resolution
			p.phi_res = 100;


			%%%%
			%Output the path properties
			output = p;
	end
	
end

function [stroke] = strokedef(tau)
    
    % Leg swing amplitude (consider a single sine gait for now)
	amp_alpha = 0.25;
    
    % Compute the leg swing shapes directly since they are controlled
    stroke = [amp_alpha*sin(tau) + 0.3, amp_alpha*cos(tau) - 0.3]; % These are your leg 
    % angles as a function of time.
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % COMPUTE SPRAWL -- not independently controlled ~~~~~~~~~~~~~~~~~~~~~~
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % We shall solve this as a differential equation where we know sprawl
    % velocity as a function of current sprawl and limb angles.

    % Define the initial sprawl
    g_0 = pi/4; % Here, we assume the robot starts from it's initial
    % mean sprawl configuration

    % Get the sprawl solution:
    [T,Y] = ode45(@sprawl_velocity, [0, 2*pi], g_0);

    % Now, get the sprawl for the times requested:
    stroke(:,end+1) = interp1(T,Y,tau,'spline'); % add another column
    
    % Check the sysf file: 
    % "\GeometricSystemPlotter\UserFiles\
    % GenericUser\Systems\sysf_sysf_hybridC_2c2r_multidisk.m"

end

function gd = sprawl_velocity(t,g)

    % Leg swing amplitude (consider a single sine gait for now)
	amp_alpha = 0.25;

    % Define the map from limb rotational velocity to sprawl velocity
    Omega_1p_2ndrow = @(a1,a2,g) [0.5*(sin(a1-a2+2*g)+cos(a1)-cos(a1+2*g)), -0.5*(sin(a1-a2+2*g)+cos(a2)-cos(a2-2*g))]*...
        (1/(sin(2*g) - sin(a1-a2+2*g) + cos(a1+2*g) + cos(a2-2*g)));

    % Get the limb swing angle and angular velocities:
    a1 = amp_alpha*sin(t) + 0.3; a2 = amp_alpha*cos(t) - 0.3;
    ad = [amp_alpha*cos(t), -amp_alpha*sin(t)]';

    % Obtain the sprawl velocity:
    gd = Omega_1p_2ndrow(a1,a2,g)*ad; % scalar (2x1)(1x2)

end
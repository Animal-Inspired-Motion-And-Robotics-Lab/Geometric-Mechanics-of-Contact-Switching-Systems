function output = shchf_diametricallyOppLegs_SwingAnglesOnly(input_mode,~)

	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end	
	
	switch input_mode
		
		case 'name'
			
			output = '2D Leg Swing Angle: Sinusoidal Alpha, Amp = pi/2';
			
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

function [stroke] = strokedef(t)

	t = -t(:)';

	Rot=sqrt(2)/2*[1 -1;1 1];
	a=pi/2;

	stroke=(Rot*[-a*cos(t);-a*sin(t)] + [0, 0]')'; % [0, pi]
    % The diametrically opposite nature of the gait is capture in the input
    % to the metric and connections functions. Check the sysf file: 
    % "C:\Users\hario\Documents\GitHub\GeometricSystemPlotter\UserFiles\
    % GenericUser\Systems\sysf_two_pin_legged_rigidbody_debugVer.m"

end
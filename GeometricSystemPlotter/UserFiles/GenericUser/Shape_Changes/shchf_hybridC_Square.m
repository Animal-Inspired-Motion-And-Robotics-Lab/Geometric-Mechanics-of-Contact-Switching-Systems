function output = shchf_hybridC_Square(input_mode,pathnames)

	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end	
	
	switch input_mode
		
		case 'name'
			
			output = 'Swing Angle and Adhesion: Square Gait, Amp = pi/4';
			
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

	amp_alpha = 1.1; % pi/4, 1.1, 1.5

    % Create a negative sign flag:
    n = @(phi) -1*(phi < 0);

    % Normalize the phi-range:
    phiN = @(phi) mod(phi,n(phi)*2*pi);

    % Create some inline functions:
    triangularAlpha = @(t) ((2*amp_alpha/pi)*(n(phiN(t)).*phiN(t) - pi)...
        - amp_alpha).*(n(phiN(t)).*phiN(t) >= pi)...
        + ((-2*amp_alpha/pi)*n(phiN(t)).*phiN(t) + amp_alpha).*(n(phiN(t)).*phiN(t) < pi);
    stepC = @(t) 2*double(n(phiN(t)).*phiN(t) < pi) -1; %double(n(phiN(t)).*phiN(t) < pi) % 2* %-1
    
    % Now, we get the function
    stroke = [triangularAlpha(t), stepC(t)];

    % The diametrically opposite nature of the gait is capture in the input
    % to the metric and connections functions. Check the sysf file: 
    % "C:\Users\hario\Documents\GitHub\GeometricSystemPlotter\UserFiles\
    % GenericUser\Systems\sysf_two_pin_legged_rigidbody_SmoothModeChange.m"

end
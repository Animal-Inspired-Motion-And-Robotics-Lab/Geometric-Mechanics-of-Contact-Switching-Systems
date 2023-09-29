function output = shchf_diametricallyOpp_SwingAdhesionCoords_sinusoidal_inPhase(input_mode,pathnames)

	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end	
	
	switch input_mode
		
		case 'name'
			
			output = 'Swing Angle and Adhesion: in-phase Contact, Sinusoidal Alpha, Amp = pi/2';
			
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

    phi_off = deg2rad(1); % for this case the swing and contact angle are 
                 % in-phase (with a small 1deg offset)
%     phi_off = deg2rad(179);
    
	amp_alpha = pi/2;
    
%     stroke = [amp_alpha*cos(t), (1+cos(t + phi_off))/2]; % OLD FORMULATION [0,1] c
    stroke = [amp_alpha*cos(t), cos(t + phi_off)]; % NEW FORMULATION [-1,1] c
%     stroke = [amp_alpha*cos(t), (1-sin(t))/2]; % fixed 270degs ph-lead (CCW gait)

    % The diametrically opposite nature of the gait is capture in the input
    % to the metric and connections functions. Check the sysf file: 
    % "C:\Users\hario\Documents\GitHub\GeometricSystemPlotter\UserFiles\
    % GenericUser\Systems\sysf_two_pin_legged_rigidbody_SmoothModeChange.m"

end
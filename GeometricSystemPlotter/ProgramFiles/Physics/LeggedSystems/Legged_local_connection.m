function [A, h, J,J_full,Omega] = Legged_local_connection(geometry,physics,shapeparams)
% Calculate the local connection for a set of curvature bases
%
% Inputs:
% geometry: structure defining system geometry
%      geometry.type: how the system geometry is defined 
%         (e.g., n-disk-legged, n-multidisk-legged, discrete-compliant)
%      geometry.function: map from shape variables to local backbone
%         deformation (e.g., joint angles)
%      geometry.length: total length of the walker
% physics: structure defining system physics
%      drag_ratio: ratio of lateral to longitudinal drag
%      drag_coefficient: drag per unit length
% cparams: value of shape variables




% Identify what kind of system is being calculated, and use this to specify how
% the local connection should be generated
switch geometry.type

    case 'n-disk-legged'

        % Now, we make sure there is only one pinned leg:
        if sum(isinf(physics.st_friction_coeff),'all') == 1

            physics_function = @pinned_leg_constraint;

        elseif sum(isinf(physics.st_friction_coeff),'all') > 1
            
            % If the number of pinned feet equal to the number of modes,
            % then initialize the 
            if sum(isinf(physics.st_friction_coeff),'all') == size(physics.st_friction_coeff,1)
                
%                 physics_function = @pin2slip_leg_constraint; % OLD
%                 DOMAIN -- c \in [0,1]
                physics_function = @pin2slip_leg_constraint_ext; % NEW 
                % DOMAIN -- c \in [-1,1]

            else
                
                error(['ERROR: Single disk legged system cant have more than 1 ' ...
                'pinned leg']);

            end
        else

            error('ERROR: Non-pinned contact cases not supported at the moment!');

        end


    case 'n-multidisk-legged'

        error('ERROR: Multidisk systems are not supported at the moment!');
        
end

% Call the physics function identified for the system
[A, h, J, J_full, Omega] = physics_function(geometry,physics,shapeparams);

end
function [A, h, J, J_full, omega] = HybridContact_connection_discrete_v2(geometry,~,jointangles)
% Calculate the local connection for a set of curvature bases
%
% Inputs:
%
%   geometry: Structure containing information about geometry of chain.
%       Fields required for this function are:
%
%       linklengths: A vector of link lengths defining a kinematic chain
%
%       baseframe (optional): Specification of which link on the body should be
%           treated as the base frame. Default behavior is to put the base
%           frame at the center of the chain, either at the midpoint of the
%           middle link, or at the joint between the two middle links and at
%           the mean of their orientation. Options for this field are
%           specified in the documentation of N_link_chain.
%
%       length (optional): Total length of the chain. If specified, the elements of
%           will be scaled such that their sum is equal to L. If this field
%          is not provided or is entered as an empty matrix, then the links
%          will not be scaled.
%
%   physics: Structure containing information about the system's physics.
%       Fields are:
% 
%       drag_ratio: Ratio of lateral to longitudinal drag
%
%       drag_coefficient: Drag per unit length for longitudinal direction
%
%   jointangles: Rotor angle 'alpha' - offsets between each paddle is
%       predetermined in the subchains structure defined in the
%       sysf_trefoil_lowRe_expt -- the contact variable is removed within
%       the branched chain formulation
%
% Outputs:
%
%   A: The "local connection" matrix, or "Locomotion Jacobian". This matrix
%       maps joint angular velocities ("shape velocities") to body
%       velocities of the system's base frame as 
%
%            g_b = -A * alphadot
%
%       with g_b the body velocity and alphadot the joint angular velocity.
%
%       (Note the negative sign in this equation; sysplotter code uses the
%       classical formulation of the geometric mechanics equations, which
%       includes this negative sign.
%
%   [h,J,J_full]: location and Jacobians for the links on the chain. These
%       are passed through from N_link_chain; see the documentation on that
%       function for more details
%
%   omega: The "Pfaffian constraint matrix" from which the local connection
%       A is calculated. This matrix is a linear map from system body and
%       shape velocities to net external forces acting on the base frame,
%       which must be zero for all achievable motions of the system
    
% Get the number of subchains and repeat the rotor angle -- attachment
    % of the subchains will add angle to each chain to maintain orientation
    % with the base link
    numC = numel(geometry.subchains);

    %%%%

    % First, get the contact states from the shape-space input:
    C = jointangles(2);

    %%%%
    % Second, get the positions of the links in the chain and their
    % Jacobians with respect to the system parameters
    if isfield(geometry,'subchains')
        [h, J, J_full] = branched_chain(geometry,jointangles);
    else
        [h, J, J_full] = N_link_chain(geometry,jointangles);
    end

    % Make an inline function that will compute the strength of the
    % connection -- higher at actual contact since 
    Cs = @(d) 2;

    % Create an empty Omega matrix to store the Pfaffian constraints:
    omega = repmat({nan(size(J_full{1})-1)},1,numel(J_full)/2);

    % Create a container to store the A matrix:
    A = zeros(size(J_full{1},1),numel(jointangles));

    % Now we need to modify the Jacobian since there is only one rotor
    % angle, so J will have 1 column and J_full will have 4 columns. Since
    % only even numbered Jacobians indicate maps from actual leg-links, we
    % only get the pfaffian constraints them.
    for idx = 1:numel(J_full)
        % Get the standard derivative of h -- J:
        J{idx} = [sum(J{idx},2),zeros(3,1)]; % The zeros correspond to contact variation
        % Now, get the contact Jacobian J_full:
        J_temp = J_full{idx}(:,4:end); J_temp = sum(J_temp,2);
        J_full{idx}(:,4) = J_temp; J_full{idx}(:,5) = zeros(3,1);
        J_full{idx} = J_full{idx}(:,1:5);
        % Check if this is a rotor link
        if mod(idx,2) == 0
            % Get the index to work with -- check if it is a 3 or two contact
            widx = idx/2;
            % Contact index (independent because you can contact <=3 at the moment)
            % Get the contact value
            % system:
            f = contactMap(C,widx,numC);
            % Get the Omega matrix for this leg's contact mode:
            omega_temp = J_full{idx}; omega_temp(:,3) = [0, 0, -1]';
            omega_temp(:,4:end) = f*omega_temp(:,4:end);
            omega{idx/2} = omega_temp;
            % Get the A matrix and zero out the last row (refer
            % pin2slip_constraint.m file for more details)
            A_temp = omega_temp(:,1:3)\omega_temp(:,4:end);
            A_temp(3,:) = 0*A_temp(3,:);
            % Add the result to the A matrix to get the interpolated
            % connection:
            A = A + A_temp;
        end
    end

    % Finally modify the strngth of the connection:
    A = Cs(C)*A;

end


%%% interpolation function - version 2 (contact connects on itself)

function fi = contactMap(c,i,n)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In V2, first and last contact states are formed by 1. All other
    % contact states are intermediate.

    % Get the spacing (V2):
    s = 2/n;

    % Get the sinusoidal term coefficient:
    coeff = pi/s;
    
    % Switch the case based on which interpolation function is needed.
    switch 1


        case i == 1 % first contact state
            
            if c <= -1
                fi = 1;
            elseif c <= -1+s && c > -1
                fi = ( 1 + cos(coeff*(c+1)) )*0.5;
            elseif c > -1+s && c <= -1 + (n-1)*s
                fi = 0;
            elseif c > -1 + (n-1)*s && c <= 1
                fi = ( 1 - cos(coeff*(c + 1 -(n-1)*s)) )*0.5;
            else % greater than or equal to 1 case
                fi = 1;
            end

        case i ~= 1 % all other contact states

            if c <= -1 + (i-2)*s
                fi = 0;
            elseif c <= -1 + i*s && c > -1 + (i-2)*s
                fi = ( 1 - cos(coeff*(c + 1 - (i-2)*s)) )*0.5;
            elseif c > -1 + i*s
                fi = 0;
            end


    end



end

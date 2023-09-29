function [A, h, J, J_full, omega] = pin2slip_leg_constraint(geometry,~,shapevars)
% Calculate the local connection for a disk legged system.
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
%           frame at the center of the chain -- in our case it is at the
%           center of the rigid body.
%
%       length (optional): Total length of the chain. If specified, the elements of
%           will be scaled such that their sum is equal to L. If this field
%          is not provided or is entered as an empty matrix, then the links
%          will not be scaled.
%
%   physics: Structure containing information about the system's physics.
%       Fields are:
% 
%               static friction coefficient: For this parameter to have a
%               finite value, we need lateral in-plane torque/force
%               estimates from the actuator to establish legs that are
%               pinned to the ground vs dragging along. Currently, this is
%               set to have an 'INF' value for the legs pinned to the
%               ground, zero for the dragging ones, and NaNs for the rigid
%               body link(s) that are not contacting the ground.
% 
%               coulomb friction coefficient: This is the constant offset
%               term in the kinetic friction acting on ground-contacting 
%               moving legs. The most straightforward way to modulate this
%               would be to incorporate some form of adhesion component to
%               the feet in the model (for instance dry adhesion
%               approximated based on feet geometry or electroadhesion,
%               etc) -- in the case of electroadhesion, we can change the
%               normal component of the reaction force and hence obtain a
%               larger friction cone.
% 
%               viscous friction coefficient: This is the most
%               straightforward of terms and we consider the opposing force
%               to be linear in the contacting leg's translational
%               velocity.
%
%   shapevars: Leg hip-angles and contact states
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
    
    % Get the number of leg pairs:
    numLegPair = numel(shapevars)/2;

    % Initialize the right-side leg angles:
    jointangles(:) = [shapevars(1:numLegPair); pi + shapevars(1:numLegPair)];

    % Get the contact variables and check if there are any contact 
    % variables with values < 0  or greater than 1, then limit them at 0,1
    % and reassign them to shapevars --  we add a small padding that
    % slightly expands the domain so that the optimization is numerically
    % viable (maybe +-5percent of the corresponding range)
    %%%%%% UNCOMMENT TWO LINES BELOW IF YOU WANT TO DO DOMAIN CONTROL:
    c = shapevars(numel(shapevars)/2 + 1:end); 
    c(c < -0.05) = -0.05; c(c > 1.05) = 1.05;
    shapevars(numel(shapevars)/2 + 1:end) = c;
    % Switch off A outside the domain condition:
%     alpha = jointangles(1:numel(jointangles)/2); % Get the first half of the joint angles
%     domFlagA = sum(alpha < -pi/2 | alpha > pi/2) == 0;
%     domFlagC = sum(c < -0.05 | c > 1/05) == 0;
%     domFlag = double(domFlagA || domFlagC); % if it is outside the domain

    dom_mid = 0.5; % mid of the contact input range

    %%%%
    % First, get the positions of the links in the chain and their
    % Jacobians with respect to the system parameters -- here, we know that
    % the system has leg-pairs that are diametrically opposite and that
    % gives us the number of 
    [h, J, J_full] = n_disk_legs(geometry,jointangles);

    %%%%%%%%
    % We here impose that one of the legs are pinned to the ground by
    % imposing 0 translational velocity to the leg in question and zero
    % rotational velocity for the system's identity element/ body frame.
    % This way, we end up with a hybrid local connection where the
    % shape-space velocities map to the body translation and pinned leg
    % rotationa velocities.

    % Assign the full jacobian as the omega matrix and swap the third
    % column with a [0, 0, -1] -- this comes from how we assign our
    % constraints:
    % 1) on LHS we have translational velocities for the pinned feet zero
    % with only the rotational velocity non-zero
    % 2) on the RHS we have the rotational velocity of the body-frame zero
    % Hence, we swap these to columns to get our Pfaffian matrix
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Initialize the A matrix:
    A = zeros(3,length(shapevars));
    % Initialize the Omega cell array:
    omega = zeros(size(J_full{1},1), 2*size(J_full{1},2)); 
    % We are gonna concatenate each pfaffian constraint next to each other.
    omega = repmat({omega}, 1, numLegPair);

    % Candidate functions for contact state interpolation (HAVE ONE PAIR OF
    % FUNCTIONS ACTIVE AT ANY TIME):
%     f1 = @(d,d_mid) 1-d;
%     f2 = @(d,d_mid) d;
%     f1 = @(d,d_mid) (1+cos(pi*d))/2;
%     f2 = @(d,d_mid) (1-cos(pi*d))/2;
    f1 = @(d,d_mid) 1 - 1./(1 + exp(-10*(d-d_mid)));
    f2 = @(d,d_mid) 1./(1 + exp(-10*(d-d_mid)));

    % Now, we obtain the full jacobian of this leg as the 'Omega' matrix
    % which is the pfaffian constraint on this system.
    for idx = 1:numLegPair
        % get the contact value:
        dom =  shapevars(numLegPair+idx);
        % the leg on the right
        temp1 = J_full{idx}; temp1(:,3) = [0, 0, -1]'; 
        temp1(:,4:end) = f1(dom,dom_mid)*temp1(:,4:end);
        % the leg on the left (increment once more to accomodate for the
        % rigid body link in the middle):
        temp2 = J_full{numLegPair+1+idx}; temp2(:,3) = [0, 0, -1]'; 
        temp2(:,4:end) = f2(dom,dom_mid)*temp2(:,4:end);
        % Concatenate the constraints side by side since we are diluting
        % the connections any way,
        omega{idx} = [temp1, temp2];
        % Get the columns of the local connection:
        tempA1 = temp1(:,1:3)\temp1(:,4:end);
        tempA2 = temp2(:,1:3)\temp2(:,4:end);
        % Get the specific column of the connection:
        A(:,idx) = sum(tempA1 + tempA2, 2);
        % All the columns corresponding to the contact derivative is
        % already 0.
    end


    %%%%%%%%
    % Finally, we make the third row of the A-matrix 0 since we don't want 
    % the rotational velocity of the leg-pin, rather we need to 0 the rotational
    % velocity on the body frame:
    A(3,:) = zeros(1,size(A,2));

    % Check if we are zeroing A outside the domain.
    if exist('domFlag','var')
        A = domFlag*A;
    end
    
end
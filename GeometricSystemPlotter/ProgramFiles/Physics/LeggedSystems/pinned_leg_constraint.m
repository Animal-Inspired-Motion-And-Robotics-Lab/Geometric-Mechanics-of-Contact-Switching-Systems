function [A, h, J, J_full, omega] = pinned_leg_constraint(geometry,physics,jointangles)
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
%   jointangles: Angles of the joints between the links, defining the
%       chain's current shape.
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


    %%%%
    % First, get the positions of the links in the chain and their
    % Jacobians with respect to the system parameters
    [h, J, J_full] = n_disk_legs(geometry,jointangles);

    %%%%%%%%
    % We here impose that one of the legs are pinned to the ground by
    % imposing 0 translational velocity to the leg in question and zero
    % rotational velocity for the system's identity element/ body frame.
    % This way, we end up with a hybrid local connection where the
    % shape-space velocities map to the body translation and pinned leg
    % rotationa velocities.
    %    
    % To calculate this linear hybrid local connection, we first check
    % which leg is adhered:
    idx = isinf(physics.st_friction_coeff);

    % Assign the full jacobian as the omega matrix and swap the third
    % column with a [0, 0, -1] -- this comes from how we assign our
    % constraints:
    % 1) on LHS we have translational velocities for the pinned feet zero
    % with only the rotational velocity non-zero
    % 2) on the RHS we have the rotational velocity of the body-frame zero
    % Hence, we swap these to columns to get our Pfaffian matrix
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Now, we obtain the full jacobian of this leg as the 'Omega' matrix
    % which is the pfaffian constraint on this system.
    omega = J_full{idx}; omega(:,3) = [0, 0, -1]';

    %%%%%%%%
    % Since we obtained the omega matrix by velocity constraints on the
    % pinned foot-tip and body-frame, below we show how to go from Omega to
    % the local connection A
    %
    % [0] = omega * [g_b; alphadot]
    %
    % For our locomotion studies, we will be specifying shape velocities and
    % mapping them to the body velocities that satisfy the constraint
    % equation. To construct this map, we split omega into the blocks that
    % act on g_b and alphadot,
    %
    % [0] = [omega_g omega_alpha] * [g_b; alphadot],
    %
    % pull the body velocity terms to the left,
    %
    % -(omega_g * g_b) = omega_alpha * alphadot,
    %
    % and then invert the omega_g term,
    %
    % g_b = - ( inv(omega_g) * omega_alpha ) * alphadot.
    %
    % Grouping the two omega terms then gives  the local connection as
    % 
    % g_b = - A * alphadot
    %
    % (Where the notational choice to leave the negative sign out of the
    % grouping is an historical artifact from the history of geometric
    % mechanics)
    
    % For this calculation, omega_g is the first three columns of omega,
    % and omega_alpha is the remaining columns
    A = omega(:,1:3)\omega(:,4:end);

    % Now, we make the third row of the A-matrix 0 since we don't want the
    % rotational velocity of the leg-pin, rather we need to rotational
    % velocity on the body frame:
    A(3,:) = zeros(1,size(A,2));
    
end
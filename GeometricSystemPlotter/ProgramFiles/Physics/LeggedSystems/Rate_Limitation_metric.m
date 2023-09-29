function Mp = Rate_Limitation_metric(geometry,physics,jointangles)
% Calculate the power dissipation metric for an n-legged link system based
% on the rate limits on the system's actuators
%
% Inputs:
%
%   geometry: Structure containing information about geometry of chain.
%       Fields required for this function are:
%
%       linklengths: A vector of link lengths defining a branched chain
%
%       baseframe (optional): Specification of which link on the body should be
%           treated as the base frame. Default behavior is to put the base
%           frame at the center of the chain, either at the midpoint of the
%           middle link which is the rigid body in our case.
%
%       length (optional): Total length of the chain. If specified, the elements of
%           will be scaled such that their sum is equal to L. If this field
%          is not provided or is entered as an empty matrix, then the links
%          will not be scaled.
%
%   physics: Structure containing information about the system's physics.
% 
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
%               Each row in the coefficient array highlights the different
%               modes we are interested in and the connection establishes
%               as a smooth transition from one of these modes to the other
%               bu introducing a connection vector
% 
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
%
% Output:
%
%   Mp: A matrix (m x m matrix, where m is the number of elements in
%       the joint angle input). This matrix encodes the resistance of the
%       viscous medium to motion of the joints. It acts as a linear map
%       from joint angular velocity to torques acting on the joints, and as
%       a quadratic map from joint angular velocity to power being
%       dissipated through the system. 
%
%       If Mp is taken as a metric tensor, pathlengths calculated from it
%       have lengths in units of "sqrt(power)*time". This pathlength is a
%       pacing-indepenent cost for the motion, and represents a "minimum
%       cost" for following the path, which is achieved for a trajectory
%       that follows the path at constant speed as measured by the metric,
%       which itself corresponds to constant power. Given two paths, the path
%       with a shorter length according to this metric can always be
%       traversed more quickly at a given power, or with a lower power in a
%       given time.


    %%%%%%%
    % 	To calculate this metric, first get the Local Connection, link
    % 	configurations, and full Jacobians for the links

    [A,~,~,~] = Legged_local_connection(geometry,physics,jointangles);    
    
    % Metric is the identity matrix acting on the shape-velocities
    % (weighted by the values given in the mDiag array:
    if isempty(physics.mDiag) % just use identity
        Mp = eye(size(A,2)); % positive definite rate limit metric tensor
    else
        if numel(physics.mDiag) ~= size(A,2)
            error('ERROR: The number of entries should match the size of the shape-space!');
        else
            Mp = diag(physics.mDiag);
        end
    end
%     Mp = eye(size(A,2));
    
end
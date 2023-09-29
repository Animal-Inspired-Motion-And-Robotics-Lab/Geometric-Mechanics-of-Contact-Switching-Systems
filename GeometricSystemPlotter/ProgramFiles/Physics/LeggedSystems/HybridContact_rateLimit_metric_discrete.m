function Mp = HybridContact_rateLimit_metric_discrete(~,physics,~)
% Calculate the power dissipation metric for an n-link chain under a
% viscous resitive-force model
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
%   shapeparams: Rotor angle 'alpha' and the contact variable 'c'
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
    % Let's just use a diagonal matrix to imply rate limits in each
    % direction:
    Mp = diag(physics.rateLim);
    
end
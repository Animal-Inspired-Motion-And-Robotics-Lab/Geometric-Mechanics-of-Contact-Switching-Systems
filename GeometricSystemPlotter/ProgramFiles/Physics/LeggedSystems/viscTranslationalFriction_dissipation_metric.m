function Mp = viscTranslationalFriction_dissipation_metric(geometry,physics,jointangles)
% Calculate the power dissipation metric for an n-legged link system based
% translational viscous force interaction at the leg tips
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
%           For this function, we only compute the viscous contributions to
%           the power dissipation of the system (isotropic)
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

    [A,~,~,J_full] = Legged_local_connection(geometry,physics,jointangles);
    
    %%%%%%%%
    % Now calculate the metric contribution to each leg-link

    % Pre-allocate storage for the metric contributions. Each contribution
    % is m x m (square matrix with as many rows/columns as there are shape
    % variables), and there is one contribution per leg-link that is not
    % pinned to the ground.
    link_metrics = repmat({zeros(numel(jointangles),numel(jointangles))},size(J_full));
    
    % Get the Metric contribution from each leg-link. This is the drag 
    % matrix acting on an individual link, pulled back onto the space of 
    % joint angles by the total Jacobian (including locomotion effects) 
    % from joint angular velocity to the motion of the link.
    for idx = 1:numel(link_metrics)
        
        % Check if this is a dragging leg before computing the dissipation
        if physics.st_friction_coeff(idx) == 0

            link_metrics{idx} = legged_transViscDrag_metric(J_full{idx},...% Jacobian from body velocity of base link and shape velocity to body velocity of this link
                              A,...                                        % Local connection (i.e. Jacobian from shape velocity to body velocity of base link)
                              physics.viscous_friction_coeff);             % Bulk drag coefficient 

        end

    end

    %%%%%%%
    % Finally, sum the metric contributions of the links to get the metric
    % for the whole system.
    Mp = sum(cat(3,link_metrics{:}),3);
    
end

function link_metric = legged_transViscDrag_metric(J_full,A,c)
% Calculate the contribution to the power metric from a link. This is the
% drag resistance on the link, pre- and post-multiplied by the Jacobian
% from shape velocity to this link's body velocity

	%%%%%%%
    % Local drag on a link.
    %   This is a a scaled by c, identity matrix transform on the
    %   translational velocities obtained from the full jacobian, and the
    %   final main diagonal entry is '0' since leg foot point-contact
    %   rotational velocities
    
    % (Note that drag terms need to be positive here to make this a
    % positive-definite matrix that we can use as a metric)
    drag_matrix = [1 0 0;
                   0 1 0;
                   0 0 0]*c;

    %%%%%%%%%%%
    % Calculate the jacobian from shape variables to body velocity of this
    % point on the system:
    
    % -A maps shape velocity to body velocity of the system, so augmenting
    % -A with an identity matrix produces a map from shape velocity to body
    % and shape velocity of the system
    J_intermediate = [-A; eye(size(A,2))];
    
    % J_full maps the system's body and shape velocity to the body velocity
    % of this link.
    J_total = J_full * J_intermediate;
    
    % Pre- and post-multiplying J_total by the drag matrix pulls the drag
    % matrix back from being a metric on the link's body velocities to
    % being a metric on the system's shape velocities
    link_metric = J_total' * drag_matrix * J_total;

    % Now, we check if the metric is of the type double or is symbolic --
    % if symbolic, we can simplify it before returning:
    if isa(link_metric, 'sym')

        link_metric = simplify(link_metric, 'Steps', 10);

    end

        %%%%%%%%%%%
    
%     % Let's now create the shape-space vel symbolic arrays based on the 
%     % number of hip-joint angles available to us (number of columns in A is
%     % the size of these)
%     alpha_dot = sym('alpha_dot%d', [m, 1], 'real');
% 
% 	
%     % Velocity of our rigid-body COM link is same as the first two rows of 
%     % the local connection, and the rotational velocity is 0. This comes
%     % from our single leg pin connection connection.
%     g_circ = transMat*-A*alpha_dot;
% 
%     % Now, we obtain the velocity at the leg-tip which is a function of
%     % joint velocities and body velocity of the system and full jacobian
%     % carries this out
%     g_circ_leg = J_full*g_circ;
% 
%     % Let's now get the power loss formulation -- this will be F.v
%     % where the 'dot' stands for the dot product and F is the dissipation
%     % force, and since F is linear in v, we have F.v being v'*D*v where D
%     % is the drag matrix we obtained earlier.
%     P = g_circ_leg'*drag_matrix*g_circ_leg;
% 
%     % Now, P will be a 2nd order polynomial expression in alpha_dot vector 
%     % and we want to extract the coefficients of various terms that will 
%     % help us artificially reconstruct the link metric
%     [coeffs_IP,  terms_IP] = coeffs(P, alpha_dot, 'All');
% 
%     % Lets switch off this warning for the upcoming operations:
%     warning('off', 'symbolic:sym:isAlways:TruthUnknown');
% 
%     % Initialize the link-metric matrix:
%     link_metric = zeros(m, m);
% 
%     % Let's cycle through each term in a nested loop and obtain our
%     % coefficients (since the metric tensor is symmetric we can obtain the
%     % coefficients for the upper-triangular matrix and then split it
%     % equally between the off-diagonal terms)
%     for idx1 = 1:m
% 
%         for idx2 = idx1:m
% 
%             % Obtain the terms we want coefficients for:
%             symComp = str2sym(['alpha_dot_',num2str(idx1)])*...
%                 str2sym(['alpha_dot_',num2str(idx2)]);
%             
%             % Check if it is a main diagonal term or not:
%             if idx1 == idx2
%                 
%                 % Obtain the coefficient:
%                 link_metric(idx1, idx2) =...
%                     sum(isAlways(terms_IP == symComp).*coeffs_IP,'All');
% 
%             else
%                 
%                 % Populate both off-diagonal terms:
%                 link_metric(idx1, idx2) =...
%                     0.5*sum(...
%                     isAlways(terms_IP == symComp).*coeffs_IP,'All');
% 
%                 link_metric(idx2, idx1) = link_metric(idx1, idx2);
% 
%             end
% 
%         end
% 
%     end
    
end
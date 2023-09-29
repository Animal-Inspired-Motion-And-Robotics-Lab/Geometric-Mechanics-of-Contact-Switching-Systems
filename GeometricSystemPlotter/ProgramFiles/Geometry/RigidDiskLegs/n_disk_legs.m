function [h, J, J_full] = n_disk_legs(geometry,shapeparams)
% Build a backbone for a legged system (a special branched chain), 
% specified as a vector of link lengths and the joint angles between them.
%
% Inputs:
%
%   geometry: Structure containing information about geometry of chain.
%       Fields are:
%
%       linklengths: A vector of link lengths defining a kinematic chain
%
%       baseframe (not-supported): Specification of which link on the body should be
%           treated as the base frame. Default behavior is to put the base
%           frame at the center of the chain, either at the midpoint of the
%           middle link, or at the joint between the two middle links and at
%           the mean of their orientation. Options for this field are:
% 
%               'centered' :    Default behavior
%               BELOW NOT SUPPORTED ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%               'tail' :        Lowest-numbered link is the base frame
%               'tail-tip' :    Start of lowest-numbered link is the base frame
%               'head' :        Highest-numbered link is the base frame
%               'head-tip' :    End of the highest-numbered link is the base frame
%               numeric :       Specify a link number to use as a base frame
%               sysf_           Pull minimum-perturbation coordinates from a
%                                       sysf_ file. Argument should be the name of
%                                       a system in the current UserFiles folder
%               'start' :       Modifier on a link specification (e.g., 
%                   {2,'start'} to put the frame at the proximal end of the
%                   specified link
%               'end' :         Modifier on a link specification (e.g., 
%                   {head,'end'} to put the frame at the distal end of the
%                   specified link
%               transform:      A 3x3 SE(2) matrix giving the position
%                   of the base frame. This can be a standalone entry or a
%                   modifier on any other frame.
%
%       modes (optional): Option to map input "jointangles" across
%           multiple links in a chain (which can have more joints than the
%           provided number of "jointangles". Should be a matrix in which
%           each column is the set of joint angles on the full chain associated
%           with a unit value of the corresponding "jointangle" input.
%
%
%       length (optional): Total length of the chain. If specified, the elements of
%           will be scaled such that their sum is equal to L. If this field
%           is not provided or is entered as an empty matrix, then the links
%           will not be scaled.
%
%   shapeparams: A vector of the angles between the links. Not necessarily
%               one shorter than the number of links.
%
%   
%
%
% Outputs:
%
%   h : Locations of the chain links relative to the selected base frame,
%           and lengths of the links
%
%       h.pos: The link locations are stored in an Nx3 array, with
%           each row a link and the columns corresponding to x,y,theta
%           values for that link. These can be converted to stacks of SE(2)
%           matrices via vec_to_mat_SE2, and back to vectors by
%           mat_to_vec_SE(2)
%   
%   	h.lengths : Vector of linklengths. This is the original link length
%           specification passed into the function, scaled by L, and stored
%           as a column of values.
%
%   J : Jacobians from joint angle velocities to velocity of each link
%           relative to  the base frame (the standard derivative of h).
%           These Jacobians are stored in a cell array, with one matrix per
%           cell.
%
%   J_full : Full Jacobians from body velocity of the base frame and joint
%           angle velocities to body velocity of each link relative to a
%           non-moving frame.     

    %%%%%%%%%%%%
    % Input parsing
    
    % If no length is specified, do not scale the links
    if ~isfield(geometry,'length') || isempty(geometry.length)
        L = sum(geometry.linklengths);
    else
        L = geometry.length;
    end
    
    % If no modes are specified, use an identity mapping for the modes
    if ~isfield(geometry,'modes') || isempty(geometry.modes)
        modes = eye(numel(shapeparams));
    else
        modes = geometry.modes;
    end
    
    % Force linklength and shapeparam vectors to be columns, and normalize 
    % link lengths for a total length of 1.
    linklengths = geometry.linklengths(:)/sum(geometry.linklengths)*L;
    shapeparams = shapeparams(:);
    
    
    % Expand jointangles from specified shape variables to actual joint 
    % angles by multiplying in the modal function
    jointangles = modes*shapeparams;
    
    %%%%%%%%%%%%
    
    % Prevent Matlab from playing tricks with imaginary numbers on symbolic
    % inputs and from complaining about assumptions on constants
    if isa(jointangles,'sym')
            
        assume(symvar(jointangles),'real');
    
        warning('off','symbolic:sym:sym:AssumptionsOnConstantsIgnored')
    
    end
    
    %%%%%%%%%%%%%%
    % Get some basic information about the chain whose kinematics we're
    % calculating
    
    % Number of links and joints
    N_links = numel(linklengths);
    M_joints = numel(jointangles);
    
    % Location of the middle link --  acts as our rigid body
    Rigid_link_id = ceil((N_links+1)/2);
    
    % Number of links that are legs in the system:
    N_legs = N_links - 1;
    
    % % Decide if there is an even or odd number of joints
    % N_odd = mod(N_links,2);
    
    % Let's do a check to see if the number of joint angles is equal to the
    % number of legs:
    if N_legs ~= M_joints
        error('ERROR: Number of joint angles should equal the number of legs');
    end
    
    %%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%
    % Core kinematics code
    
    % Introduce a transforms structure:
    h = [];
    
    % Firstly, the link lengths can directly be ported to the introduced 
    % struct
    h.lengths = linklengths;
    
    % Create empty jacobian cell-array-matrix hybrids:
    J = zeros(3, M_joints); J = repmat({J}, 1, N_links);
    J_full = zeros(3, M_joints + 3); J_full = repmat({J_full}, 1, N_links);
    
    % The center of the rigid body is the body frame -- identity element.
    % Let's find the translations to the legs from this link. The proximal 
    % ends of the leg-links are attached to the center of the body-frame. 
    % 
    % Each leg is only affected by the corresponding joint angle.
    % 
    % Hence, we have a straightforward transform to the distal-end of the 
    % legs as defined in the anonymous function below,
    leg_distal_transform = @(a,l) l*[cos(a), sin(a), a];
    
    % Iterate over the legs to obtain the Transforms and the Jacobians
    for idx1 = 1:N_links
    
        % Set the Jacobian compute flag to false:
        Jf_flag= false;
        
        % Check if it is the body or the leg and compute the transform
        if idx1 == Rigid_link_id
    
            leg_flag = false;
            
            % If this is not a leg, compute the transform now:
            h.pos(idx1,:) = leg_distal_transform(0,0);
    
            % The first 3x3 matrix of the full jacobian:
            J_full{idx1}(:,1:3) = Adjinv(h.pos(idx1,:));
    
            % Set the jacobian compute flag to true:
            Jf_flag = true;
    
        else 
    
            leg_flag = true;

            % The transform if it is a leg with the corres. joint angle:~~~
            % Check which is the corresponding angle:
            if idx1 < Rigid_link_id
                idx2 = idx1;
            else
                idx2 = idx1 - 1;
            end
            % Do the compute
            h.pos(idx1,:) =...
                leg_distal_transform(jointangles(idx2),linklengths(idx1));
    
        end
    
        % Iterate over each hip-joint
        for idx2 = 1:M_joints
    
            % Check if this is the corresponding leg hip-joint:
            if idx1 > Rigid_link_id
                joint_flag = (idx2 == idx1 - 1);
            elseif idx1 < Rigid_link_id
                joint_flag = (idx2 == idx1);
            else % if it is the rigid body link, joint flag is always 0
                joint_flag = false;
            end
    
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
            % Obtain the standard derivate of h - Jacobian:~~~~~~~~~~~~~~~~
    
            J{idx1}(:, idx2) = Adjinv(h.pos(idx1,:))*...
                ([0, 0, 1]*double(joint_flag)*double(leg_flag))';
    
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
            % Compute the full jacobian:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
            if ~Jf_flag
    
                % The body-vel component of the full jacobian and set the 
                % flag:
                J_full{idx1}(:, 1:3) = Adjinv(h.pos(idx1,:));
                Jf_flag = true;
    
                % Insert the joint jacobian:
                J_full{idx1}(:, 3 + idx2) = J{idx1}(:, idx2);
    
            else % If the body vel map is already computed
                
                J_full{idx1}(:, 3 + idx2) = J{idx1}(:, idx2);
    
            end
    
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end
    
    end
    
    % If we did have symbolic inputs, let's simplify the Jacobians
    if isa(jointangles,'sym')
            
        J = simplify(J, 'Steps', 10); 
        J_full = simplify(J_full, 'Steps', 10);
    
    end

end
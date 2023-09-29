function traj = ComputeSingleDisk(sim,v)
    %%%%% This function animates a single-disk, contact-swtching robot, and
    %%%%% returns the generated trajectory.
    %%%%% Also, the legs are coupled and this function only supports ONE
    %%%%% CONTINUOUS SHAPE -- hence, one rotor to control the robot, and
    %%%%% the legs are one unit long.
%     Inputs: sim, v -- both structures
%           sim: includes the robot parameters for animation
%           v: includes video parameters such frame-rate, video name, etc.

%     % Check if the video structure is provided
%     if nargin < 2
%         r = false; % Set the recording flag to false
%     else
%         r = true; % if provided, set the recording flag to true
%     end
    
    % Unpack the structures --  simulation
    phi_pin_ord = sim.phi_pin_ord; % The order in which the legs are pinned
    % can be arbitrarily long since revisiting manifolds is allowed.
    
    % If the first and last submanifolds are the same, then no switching
    % happens at the end, and that changes the time period of the gait
    if phi_pin_ord(1) == phi_pin_ord(end)
        midF = true;
    else
        midF = false;
    end

    phi_range = sim.phi_range; % 'phi_range' 2x. vector where the first row
                               % corresponds to the lower value and the
                               % second row corresponds to the higher
                               % value. arbitrary columns based on number
                               % of subaminfold switches within a gait
                               % cycle.

    T = sim.T; % 'T' encodes the time-period of system into a 2x1
               % array with the first entry being swinging time and the
               % second being switching time

    off = sim.off; % 'off' encodes the relative angular offsets
                   % between the legs encoded by mechanical coupling

%     phi_c = sim.phi_c; % color of the trajectory

    phi_tf= sim.phi_tf;  % the final time for the animation to run upto
    % Make sure this is greater than  or equal to atleast one time period

    % Unpack the structures --  video
%     v_name = v.v_name; % name of the video
    frame_r = v.frame_r; % number of frames per second in the video
    % NOTE that this input needs to be high enough

    % Get the number of legs
    numL = numel(off);
    % If you want to be more careful, check if the num leg dimensions in
    % each data is the same and generate an error if not.

    % Create some time variables useful for later
    T_leg = T(1) + T(2); % To complete the shape changes associated with one leg
    if midF
        T_phi = T_leg*(length(phi_pin_ord)-1) + T(1); % Swing + Switch time for each leg
                                       % EXCEPT THE LAST SWITCH
    else
        T_phi = T_leg*length(phi_pin_ord); % Swing + Switch time for each leg
    end
%     T_alpha = T(1); T_beta = T(2);
    if phi_tf < T_phi
        error(['ERROR: The final time value should be atleast as big as the' ...
            'provided time period.']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% COMPUTE TRAJECTORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get the time vector -- based on the frame-rate and final time
    t = linspace(0,phi_tf,floor(phi_tf*frame_r));

    % Create a gait-phase vector that we shall to compute the net
    % displacement:
    t_tau = t(t <= T_phi); % strictly less than since at t == T_phi it is
                          % the next gait cycle.

    % Initialize the contact trajectory -- phi_beta
    phi_beta = nan(size(t_tau));
    % Taking a slightly different approach from the paper, here beta takes
    % the values corresponding to the leg's index (easier to implement this
    % way)

    % Initialize the rotor trajectory -- phi_alpha
    phi_alpha = phi_beta; % same nans

%     % Rotor angle '0' position poiting vector (unit length)
%     Alpha_r = [1,0]';

    % Initialize the local connection function:
    syms a o real
    x_dot = -sin(a + o);
    y_dot = cos(a + o);

    % Initialize the COM trajectory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Xcm = nan(1,numel(t)); 
    Ycm = nan(1,numel(t));

    % Initialize the Leg positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Xpos = nan(numL,numel(t));
    Ypos = nan(numL,numel(t));
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Compute for each leg ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for i = 1:numL
        
        % Check the permute to see if there is any time spent on this
        % manifold, if yes, then proceed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if ~isempty(phi_pin_ord(phi_pin_ord == i)) %%%%%%%%%%%%%%%%%%%%%%%%

            % Get the time block(s) spent on this submanifold
            i_time_block = find(phi_pin_ord == i);

            % iterate over each submanifold occurence of this specific leg
            for j = 1:length(i_time_block)
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                % Compute the contact trajectory -- beta ~~~~~~~~~~~~~~~~~~
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                % Find the next submanifold
                if phi_pin_ord(i_time_block(j)) == phi_pin_ord(end) && ~midF...
                        && j == length(i_time_block) % make sure last
                    % if last submanifold block and new submanifold in the
                    % next gait cycle
                    nextS = phi_pin_ord(1);
                elseif phi_pin_ord(i_time_block(j)) == phi_pin_ord(end) && midF...
                        && j == length(i_time_block)
                    % if last submanifold block and continues next gait
                    % cycle
                    nextS = nan;
                elseif phi_pin_ord(i_time_block(j)) ~= phi_pin_ord(end) || j < length(i_time_block)
                    % if not the last submanifold block
                    nextS = phi_pin_ord(i_time_block(j) + 1);
                    if nextS == i % error out if it is the same as this submanifold
                        error(['ERROR: Make sure subamniolds arent stacked in ' ...
                            'phi_pin_ord!']);
                    end
                end
                
                % Before the end of the it's swinging time and after the 
                % last switch; then linspace during the switch to the next
                % submanifold.
                i_swing = (t_tau < (i_time_block(j)-1)*T_leg + T(1)) &...
                    (t_tau >= (i_time_block(j)-1)*T_leg);

                % Check if you need a switching phase
                if ~isnan(nextS)

                    % Compute the switching indices
                    i_switch = (t_tau < i_time_block(j)*T_leg ) &...
                    (t_tau >= (i_time_block(j)-1)*T_leg + T(1));

                    % Assign the switching in shape (padded after)
                    temp = linspace(i,nextS,sum(i_switch)+1);
                    phi_beta(i_switch) = temp(1:end-1); % padding removed

                end
                

                % Compute the time indices when SWINGING and SWITCHING 
                % happens -- needed later for computing displacement 
                % (switching too for computing leg positions)
                i_swing_time = find(i_swing == 1);
                
                % Assign the current submanifold
                phi_beta(i_swing) = i;
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                % Compute the trajectory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                % Get the velocities for the current leg
                x_dot_curr = subs(x_dot,o,off(i));
                y_dot_curr = subs(y_dot,o,off(i));

                % Get the displacements from previous submanifold
                % shape-changes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ord_pre = phi_pin_ord(1:i_time_block(j)-1);
                % Initialize an empty displacement containers
                delX_pre = 0; delY_pre = 0;
                if ~isempty(ord_pre) % if not empty compute the 
                                     % displacements
                    
                    % Iterate over each term and add it to integral
                    for k = 1:length(ord_pre)
                        delX_pre = delX_pre...
                            -double(int(subs(x_dot,o,off(ord_pre(k))), a,...
                            [phi_range(1,k) phi_range(2,k)]));
                        delY_pre = delY_pre...
                            -double(int(subs(y_dot,o,off(ord_pre(k))), a,...
                            [phi_range(1,k) phi_range(2,k)]));
                    end

                end

                % Add this up to the trajectory containers:
                Xcm(i_swing) = delX_pre; Ycm(i_swing) = delY_pre;
                
                % Get the current displacement ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                % Create a vector for rotor angles for this time
                alpha_curr = interp1(linspace(0,T(1)),...
                    linspace(phi_range(1,i_time_block(j)),...
                    phi_range(2,i_time_block(j))),...
                    t_tau(i_swing) - (i_time_block(j)-1)*T_leg);
                alpha_curr = alpha_curr(:)'; % make sure it is a row vector

                % Iterate at swing instant:
                for l = 1:length(i_swing_time)
                    
                    % Compute COM traj (old disp + current swing disp)
                    Xcm(i_swing_time(l)) = Xcm(i_swing_time(l)) +...
                        -double(int(x_dot_curr, a,...
                        [phi_range(1,i_time_block(j)) alpha_curr(l)]));
                    Ycm(i_swing_time(l)) = Ycm(i_swing_time(l)) +...
                        -double(int(y_dot_curr, a,...
                        [phi_range(1,i_time_block(j)) alpha_curr(l)]));

                end
                
                % If a switching phase exists, compute the net trajectory
                % during swing and assign it
                if ~isnan(nextS)

                    % Trajectory during switching phase
                    Xcm(i_switch) = delX_pre...
                        -double(int(subs(x_dot,o,off(i)), a,...
                        [phi_range(1,i_time_block(j)) phi_range(2,i_time_block(j))]));
                    Ycm(i_switch) = delY_pre...
                        -double(int(subs(y_dot,o,off(i)), a,...
                        [phi_range(1,i_time_block(j)) phi_range(2,i_time_block(j))]));

                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % Compute the swing trajectory -- alpha ~~~~~~~~~~~~~~~
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    phi_alpha(i_swing | i_switch) = [alpha_curr,...
                        phi_range(2,i_time_block(j))*ones(1,sum(i_switch))];

                else

                    % No need to compute trajectory since it is already
                    % computed.
                    
                    phi_alpha(i_swing) = alpha_curr; % only swing phase
                    
                end
                
            end
            
        end
        
    end
    
    % Compute the net displacement over a cycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    netX = 0;  netY = 0;
    for i = 1:length(phi_pin_ord) % each chronological submanifold block
        netX = netX + -double(int(subs(x_dot,o,off(phi_pin_ord(i))),a,...
            [phi_range(1,i) phi_range(2,i)]));
        netY = netY + -double(int(subs(y_dot,o,off(phi_pin_ord(i))),a,...
            [phi_range(1,i) phi_range(2,i)]));
    end

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % INTERPOLATE THE TRAJECTORY FOR THE NEXT FEW GAIT CYCLES ~~~~~~~~~~~~~
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % This method would be exposed if the frame-rate is low.

    % Now that we have the COM trajectories over the gait-cycle, let's now
    % compute net displacement over the gait cycle for final time phi_tf >
    % T_phi
    if abs(floor(phi_tf/T_phi) - phi_tf/T_phi) < 1e-5 % if the final time
                                                      % is multiple of T
        numC = floor(phi_tf/T_phi);
    else
        numC = floor(phi_tf/T_phi) + 1;  % number of gait-cycles -- 
                                     % full and 1 last partial
    end

    % Initialize empty containers for storing the shape-changes:
    Phi_alpha = nan(1,length(t)); Phi_beta = nan(1,length(t));

    % Iterate and approximate the trajectories for each subsequent cycle
    % starting from 2
    for i = 2:numC

        % Compute this only once
        if i == 2
            % The indices corresponding to the lasy cycle:
            itau_pre = (t >= (i-2)*T_phi & t < (i-1)*T_phi);
        end
        
        % Get the indices corresponding to the ith gait cycle
        itau = (t >= (i-1)*T_phi & t < i*T_phi);

        % Get the gait-phase values for this cycle
        tau = t(itau) - (i-1)*T_phi;

        % Interpolate the COM trajectories:
        intp_delX = interp1(t(itau_pre),Xcm(itau_pre),tau); %,'spline'
        intp_delY = interp1(t(itau_pre),Ycm(itau_pre),tau); %,'spline'

        % Compute the net displacement till the last cycle:
        delX_pre = (i-1)*netX; delY_pre = (i-1)*netY;

        % Finally, the trajectory for this cycle
        Xcm(itau) = delX_pre + intp_delX;
        Ycm(itau) = delY_pre + intp_delY;

        % Stack the shape changes:
        Phi_alpha(itau) = interp1(t(itau_pre),phi_alpha,tau);
        Phi_beta(itau) = interp1(t(itau_pre),phi_beta,tau);
        if i == 2 % Assign the shape-changes from the first cycle.
            Phi_alpha(itau_pre) = phi_alpha;
            Phi_beta(itau_pre) = phi_beta;
        end

    end

    if numC < 2 % Compute the trajectories separately -- if it is only one gait cycle
        Phi_alpha = phi_alpha;
        Phi_beta = phi_beta;
    end

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % COMPUTE THE LEG POSITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Iterate over the each leg
    for i = 1:numL
        
        % Compute the Leg positions -- add the rotor angle 'phi_alpha'
        % + mechanical offset 'off': 
        Xpos(i,:) = Xcm + cos(Phi_alpha + off(i));
        Ypos(i,:) = Ycm + sin(Phi_alpha + off(i));
        
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % RETURN THE TRAJECTORY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Stack vertically
    traj = [Phi_alpha; Phi_beta; Xcm; Ycm; Xpos; Ypos];

end
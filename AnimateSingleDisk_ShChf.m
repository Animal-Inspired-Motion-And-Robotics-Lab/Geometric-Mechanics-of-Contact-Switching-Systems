function AnimateSingleDisk_ShChf(traj,~,leg_c,v,off) %traj_c

    % Check if the colors for each submanifold and the offset vectors are
    % of same length -- move this check to the start or remove it if you
    % want
    if size(leg_c,1) ~= numel(off)
        error('ERROR! Offset and color arrays indicate different leg num');
    end

    % Check if the video structure is provided
    if nargin < 4
        r = false; % Set the recording flag to false
    else
        r = true; % if provided, set the recording flag to true
    end

    % Compute the number of legs from the length of the offset vector
    numL = numel(off);

    % Isolate the time vector:
    t = traj(end,:);

    % Arbitrary plot parameters: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     traj_lw = 0.5; % shape change linewidth
    
    circS = 300; % size of the current point
    fracS = 0.25; % change in size of the scatter point for non submanifold regions

    fontS_title = 40; % font size for the titles,
    fontS_labels = 40; % labels, and
    fontS_txt = 25; % text
    
    num_arrw = 15; % number of arrows per submanifold
    arrw_lw = 2; % arrow line width (TRY TO MAINTAIN ratio with arrw_scale)
    arrw_scale = 0.25; % arrow autoscale factor (0.9 default)
    arrw_bound = 20; % percentage of the range to plot the vector field over

    x_aspect = 1;
    y_aspect = 1; % aspect ratio parameters


    % Get the shape change trajectories:
    Phi_alpha = traj(1,:);
    Phi_beta = traj(2,:);

    % Unpack the video parameters: 
    frame_r = v.frame_r;
    v_name = v.shch_name;

    % Manually set the symmetric hip limits:
    ank = 90; % enter in degs

    % Define percentage for plot limits in the input directions
    pct_x = 20;
    pct_y = 10; % also moves the text around

    % Get the beta and alpha evaluation values
    beta = linspace(1,numL,numL); % using intgers to denote each submanifold
    % deviates from the contact switching work a little bit
%     alpha = linspace(-deg2rad(ank),deg2rad(ank),num_arrw);
    alpha = linspace(-(1+arrw_bound/100)*deg2rad(ank),...
        (1+arrw_bound/100)*deg2rad(ank),num_arrw); % an adapative range
        % to plot the vector field over -- fixes arrow going over the
        % range
    [alpha,beta] = meshgrid(alpha,beta);

    % Generate the single disk pin connection
    Ax = @(a,b,offset) -sin(a + offset);
    Ay = @(a,b,offset) +cos(a + offset);

    % Video stuff START
    if r
        % Create the video object
        video = VideoWriter(strcat(v_name,'.mp4'),'MPEG-4'); % creates mp4
        % Set the framerate
        video.FrameRate = frame_r;
        % Set the quality
        video.Quality = 100;
        % Open the video
        open(video);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% ANIMATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot the results in a loop
    figure('units','pixels','position',[0 0 1920 1080],'Color','w') % open the figure
    % Also, set the background to white.

    % Use S_n as to denote the submanifold levels along the discrete shape
    subM_str = '$S_' + string(linspace(1,numL,numL)) + '$';
    
    % For each translation direction x and y
    for k = 1:2

        % Call the corresponding sub-block
%         subplot(1,2,k); % next to each other
        subplot(2,1,k); % stacked mode
        
        % REMOVE DOTTED LINES- NOT WORKING OUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         % Plot the entire gait trajectory with hashed lines:
%         plot(Phi_alpha,Phi_beta,'LineWidth',traj_lw,...
%             'Color',traj_c,'LineStyle','--'); %
% 
%         % hold the plots if it is the first one.
%         hold on;

        % Iterate over each foot -- corresponding to each level 1
        % submanifold when in contact with the locomotion plane
        for j = 1:numL
            A_x = alpha(beta == beta(j));
            A_y = beta(beta == beta(j));
            if k == 1
                A_u = Ax(A_x,A_y,off(j));
            else
                A_u = Ay(A_x,A_y,off(j));
            end
            A_v = zeros(size(A_x));
            % The v-component of the vector field (component along y) 
            % doesn't change since there are no kinematics in the contact
            % switching direction -- the connection vector field is
            % singular in beta

            % Plot the vector field for the corresponding manifold
            quiver(A_x,A_y,-A_u,-A_v,'Color',leg_c(j,:),'LineWidth',arrw_lw,...
                'AutoScaleFactor',arrw_scale);
            % we have a negative term here since the connection vector
            % fields are defined as a negative map for traditionality
            % purposes - \xi = -A \dot{\alpha}
            if j == 1
                hold on;
            end

            % Add the text next to each submanifold
            if j ~= numL
                text(-(1 + (pct_x-5)/100)*deg2rad(ank),j+((pct_y+5)/100),1,subM_str(j),...
                    'fontsize',fontS_txt,'Interpreter','latex',...
                    'Color',leg_c(j,:));
            else
                text(-(1 + (pct_x-5)/100)*deg2rad(ank),j-((pct_y+5)/100),1,subM_str(j),...
                    'fontsize',fontS_txt,'Interpreter','latex',...
                    'Color',leg_c(j,:));
            end

        end

        if k == 1
            % Add labels:
            xlabel('$\alpha$','FontSize',fontS_labels, 'Interpreter','latex');
            ylabel('$\beta$','FontSize',fontS_labels, 'Interpreter','latex',...
                'Rotation',0,'VerticalAlignment','middle',...
                'HorizontalAlignment','right'); 
        end

    end
    % These vector fields DO NOT CHANGE WITH TIME and hence can be plotted
    % at the beginning.

    % Iterate over each time value
    for i = 1:numel(t)

        % Plot the current point -- same point in x and y connection dirn
        for k = 1:2

            % Call the corresponding sub-block
%             subplot(1,2,k);
            subplot(2,1,k); 
            
%             text('position',[0.5 0.5],'Interpreter','latex','String','{\boldmath$\alpha$}')
            % Get the point with which we shall check if this current point
            % belongs to the same manifold
            if i == numel(t)
                chk_pt = i-1;
            else
                chk_pt = i+1; % default if not the last point
            end
            
            % Plot the current shape as a circle and set the title
            if k == 1
                if Phi_beta(i) ~= floor(Phi_beta(i))...
                        || Phi_beta(chk_pt) ~= Phi_beta(i) % if not an int
                                                           % or is
                                                           % changing
                    p1 = scatter(Phi_alpha(i),Phi_beta(i),fracS*circS,'k','filled');
                else % if a discrete value
                    p1 = scatter(Phi_alpha(i),Phi_beta(i),circS,leg_c(Phi_beta(i),:),'filled');
                end
                title('{\boldmath$A$}$_x$','FontSize',fontS_title,...
                    'Interpreter','latex');
            else
                if Phi_beta(i) ~= floor(Phi_beta(i)) % if not an integer
                    p2 = scatter(Phi_alpha(i),Phi_beta(i),fracS*circS,'k','filled');
                else % if a discrete value
                    p2 = scatter(Phi_alpha(i),Phi_beta(i),circS,leg_c(Phi_beta(i),:),'filled');
                end
                title('{\boldmath$A$}$_y$','FontSize',fontS_title,...
                    'Interpreter','latex');
            end

            % Set the axis limits and preferences here here:
%             axis tight; % the arrows almost touch the plotbox bounds
%             axis equal tight; % still same issues as above
%             axis equal; % lots of white space -- kinda works with the limits below
%             axis equal square;
%             axis equal square tight;
            xlim([-(1 + pct_x/100)*deg2rad(ank), (1 + pct_x/100)*deg2rad(ank)]);
            ylim([1 - pct_y/100, numL + pct_y/100]);
            pbaspect([x_aspect,y_aspect,1]);

            % Set the ticks here:
            set(gca,'TickLabelInterpreter','latex');

            yticklabels({}); % use text to show the submanifolds\
%             mulAlpha = floor(4*deg2rad(ank)/pi);
%             xticks(linspace(-mulAlpha*(pi/4),mulAlpha*(pi/4),2*mulAlpha+1));
            xticklabels({});

            % Output the plot:
            drawnow();

        end
        
        if r
            % Write to video:
            writeVideo(video,getframe(gcf));
        end

        % If it is not the final frame, clear it:
        if i ~= numel(t) % this command has to be at the end.
            delete(p1); % Just delete the moving point in each axis
            delete(p2);
        end

    end
    
    % Video stuff END
    if r
        % Close the video
        close(video);
    end

end

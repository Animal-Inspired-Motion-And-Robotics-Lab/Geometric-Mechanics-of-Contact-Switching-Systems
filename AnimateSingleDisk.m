function AnimateSingleDisk(traj,traj_c,leg_c,v,limX,limY,circS,tS,off)

    % Check if the video structure is provided
    if nargin < 4
        r = false; % Set the recording flag to false
    else
        r = true; % if provided, set the recording flag to true
    end
    
    % Check if the time structure is provided:
    if nargin < 8
        timeF = false;
    else
        timeF = true;
    end

    % Unpack the time vector
    if timeF
        T_leg = tS.T_leg;
        phi_pin_ord = tS.phi_pin_ord;
    end
    T_phi = sum(T_leg);

    % Compute the number of legs:
    numL = 0.5*(size(traj,1)-5); % -4 earlier, but the time vector is added
                                % obtain just the leg position rows

    % Isolate the time vector:
    t = traj(end,:);

    % Unpack the trajectory information and video data: 
    Phi_beta = traj(2,:);
    Xcm = traj(3,:); Ycm = traj(4,:);
    Xpos = traj(5:4+numL,:); Ypos = traj(4+numL+1:4+2*numL,:);
    frame_r = v.frame_r;
    v_name = v.name;

    % Get the points where the data is nan and interpolate (TEMP FIX):
    Xcm_inan = isnan(Xcm); Ycm_inan = isnan(Ycm);
    Xpos_inan = isnan(Xpos); Ypos_inan = isnan(Ypos);
    
    Xcm(Xcm_inan) = interp1(1:length(Xcm(~Xcm_inan)),Xcm(~Xcm_inan),find(Xcm_inan));
    Ycm(Ycm_inan) = interp1(1:length(Ycm(~Ycm_inan)),Xcm(~Ycm_inan),find(Ycm_inan));
    Xpos(Xpos_inan) = interp1(1:length(Xpos(~Xpos_inan)),Xpos(~Xpos_inan),find(Xpos_inan));
    Ypos(Ypos_inan) = interp1(1:length(Ypos(~Ypos_inan)),Xpos(~Ypos_inan),find(Ypos_inan));

    % Compute the maximum dimension spread  in either x or y and use that
    % to scale the parameters.
    scaleDraw = max([diff(limX),diff(limY)])/23.5; % arbitrarily chosen value
%     scaleDraw = 1; % default scaling.

    % Get the plot parameters (change it here for now)
    lW = 1.2/scaleDraw; lW_r = 0.5/scaleDraw;
    lW_O = 0.5; alp_O = 0.2; lW_traj = 2.4/scaleDraw;
    circS = circS/scaleDraw;

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

    for i = 1:size(Xcm,2)
    
        % Plot the origin with thin gray lines
        xline(0,'k','LineWidth',lW_O,'Alpha',alp_O); 
        yline(0,'k','LineWidth',lW_O,'Alpha',alp_O);
        axis equal; % axis square equal
%         axis equal;
        hold on; set(gca, 'visible', 'off');

        % Iterate over each leg and plot the leg positions
        for j = 1:numL

            % Color container -- empty
            temp_c = zeros(numL,3);

            % Get the corresponding leg color:
            temp_c(j,:) = leg_c(j,:);
    
            % Plot the legs 
            plot([Xcm(i) Xpos(j,i)], [Ycm(i) Ypos(j,i)], 'Color', temp_c(j,:), 'LineWidth', lW);

            % Get the forward difference dervative of the contact state for
            % plotting the active contact
            Phi_beta_dot = [diff(Phi_beta) Phi_beta(1)-Phi_beta(end)];
            % since it is cyclic, we can cycle back to the first step
            % contact value.
    
            % Scatter circle to the leg tips
            if floor(Phi_beta(i)) == Phi_beta(i) && Phi_beta_dot(i) == 0
                scatter(Xpos(j,i), Ypos(j,i), circS,...
                interp1(linspace(1,numL,numL),temp_c,Phi_beta(i)),'filled');
            else
                scatter(Xpos(j,i), Ypos(j,i), circS,zeros(1,3),'filled');
            end
        
        end

        % Compute and plot the trajectory:
        if timeF
            
            % Check which gait cycle we are in:
            if t(i)/T_phi == ceil(t(i)/T_phi)
                gaitN = ceil(t(i)/T_phi)+1;
            else
                gaitN = ceil(t(i)/T_phi);
            end
            

            % Set the trajectory plotting complete flag to false
            comF = false;

            % Iterate over the gait cycles and plot the trajectories:
            for j = 1:gaitN
                
                % Get the time before this gait cycle:
                t_pre = (j-1)*T_phi;

                % Iterate over each submanifold:
                for k = 1:length(phi_pin_ord)
                    
                    % Time at the start of this submanifold.
                    if k >=2
                        t_minus = t_pre + sum(T_leg(1:k-1));
                    else
                        t_minus = t_pre;
                    end

                    % Time at the end of the manifold or current time
                    % (whichever happens earlier)
                    t_plus = t_pre + sum(T_leg(1:k));
                    if t(i) < t_plus
                        t_plus = t(i); % we plot the trajectory till this iteration
                        comF = true; 
                    end

                    % Obtain the first and last index for the trajectory
                    % plot:
                    t_idx_L = find((t >= t_minus) & (t < t_plus),1,'first');
                    if comF
                        t_idx_H = find((t >= t_minus) & (t <= t_plus),1,'last');
                    else
                        t_idx_H = find((t >= t_minus) & (t < t_plus),1,'last')+1;
                        % we do this last step to make sure these 
                        % trajectory segments are connected.
                    end

                    % Obtain the color of the trajectory plot:
                    col_temp = leg_c(phi_pin_ord(k),:);

                    % Make the trajectory plot for this section:
                    Xcm_temp = Xcm(t_idx_L:t_idx_H); 
                    Ycm_temp = Ycm(t_idx_L:t_idx_H);
                    plot(Xcm_temp, Ycm_temp, 'Color', col_temp, 'LineWidth', lW_traj);

                    % Check if the trajectory is complete and if yes, break
                    % out of the loop
                    if comF
                        break
                    end

                end
                
            end

        else % simple case of single 'k' trajectory
            Xcm_temp = Xcm(1:i); Ycm_temp = Ycm(1:i);
            plot(Xcm_temp, Ycm_temp, 'Color', traj_c, 'LineWidth', lW_traj);
        end

        % Plot a rectangle around the COM to denote the trajectory (aspect
        % ratio is fixed as of now)
        Xrect = [Xcm(i)+0.17,Xcm(i)-0.17,...
            Xcm(i)-0.17,Xcm(i)+0.17,Xcm(i)+0.17]; 
        Yrect = [Ycm(i)+0.33,Ycm(i)+0.33,...
            Ycm(i)-0.33,Ycm(i)-0.33,Ycm(i)+0.33]; 
        plot(Xrect, Yrect,'k','LineWidth',lW); % same as leg thickness

        % Plot the zero rotor angle position
        plot([Xcm(i) Xcm(i)+1], [Ycm(i) Ycm(i)], 'k--', 'LineWidth', lW_r);
    
        % Set the figure limits:
        xlim(limX); ylim(limY);
        
        % Output the plot:
        drawnow();
        
        if r
            % Write to video:
            writeVideo(video,getframe(gcf));
        end

        % If it is not the final frame, clear it:
        if i ~= size(Xcm,2) % this command has to be at the end.
            clf;
        end

    end
    
    % Video stuff END
    if r
        % Close the video
        close(video);
    end

    if nargin < 9 % the shape changes cant be plotted
        error('ERROR: you need to provide the leg offset as the 9th input');
    else
        AnimateSingleDisk_ShChf(traj,traj_c,leg_c,v,off);
    end

end
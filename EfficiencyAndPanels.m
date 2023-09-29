% This script runs an algorithm that returns an efficient SYMMETRIC gait 
% for a 2-contact system

clear all; close all; clc;

%% COLORS and FONT SIZE

% fSiz = 25;

discN = 10000;

numColors = 251;

%79.6875,49.8015,31.7156;
%37,37,37
% Custom Colormap
CUB = flipud([159.3750,99.6030,63.4312;
          239.0625,149.4045,95.1469;
          255.0000,199.2060,126.8625;
            204,204,204;
            150,150,150;
            99,99,99])*(1/255);

CUB = interp1(linspace(0,discN,size(CUB,1)), CUB, linspace(0,discN,numColors));

circ1 = [215,25,28]*(1/255); 
% circ2 = [77,172,38]*(1/255);
circ3 = [44,123,182]*(1/255);

circ2 = [96,0,220]*(1/255);

CUB_sequential = copper(numColors);

%% Plot Color and discretization:

% 2D grid discretization in each direction:
discN = 10000;

% % Custom Colormap
% CUB = flipud([158,202,225;
%             107,174,214;
%             66,146,198;
%             33,113,181;
%             8,81,156;
%             8,48,107])*(1/255);
% 
% CUB = interp1(linspace(0,discN,size(CUB,1)), CUB, linspace(0,discN,251));

% Label font size:
fontS = 10;

% Ankle limit amplitude for the rotor:
ank = pi/2 - pi/180; % keep a -pi/2 to pi/2 limit (with a 1 deg inner margin)

% Number of points to do this optimization over:
num = 100;

% COST DEFINITION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Contact switching cost:
C = linspace(0.1,1,num);
% For the 2 contact case, this cost is the time taken to disengage the
% contacting solenoid and engage the non-contacting one. Assume the units
% are in seconds.

% Swing cost -- this is the speed with which the servo motor is turning.
% Units are in rad/s -- we shall sweep a much larger range. The max
% switching time is 1s, so lets say for the max range of 2*ank rad, we need
% a max swinging time of 1s, and the slowest would be one tenth of that.
A = linspace(0.1*pi,1*pi,num); % in rad/s

% EFFICIENCY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% We assume that the gait is centered about the 0-swing position to 
% generate only net y-displacement and hence, we only need the swing 
% amplitude to get the Efficiency. Efficiency is made to take negative 
% values to suite a fmincon approach.
EffY = @(a_amp,a_cost,c_cost)...
    -2*sin(a_amp)./(c_cost + 2*a_amp/a_cost);

% Mesh the cost grid:
[a,c] = meshgrid(A,C);

% Generate the container that will store the optimal swing amplitude:
ampOptY = nan(size(a));

% Switch off the fmincon optimization results:
opts = optimoptions('fmincon','Display','none');

% Lower bound on the gait amplitude:
lb = 0.05; % in radians

% Create an array of amplitude vec:
ampV = linspace(deg2rad(1),ank,num); % we have to start with 1 deg

% Get the efficiency trajectory for each cost:
effTraj = nan([size(a),length(ampV)]);

% Let's get the optimal efficiency for each point:
EffOptY =  nan(size(a));

% Let's now find a swing amplitude that gives the y-optimal displacement:
for i = 1:length(A) % sweep over the x axis
    parfor j = 1:length(C) % % sweep over the y axis
        % Get the optimal efficiency for the given costs
        [ampOptY(i,j),EffOptY(i,j)] = fmincon(@(x) EffY(x,a(i,j),c(i,j)), 0, [], [], [],...
            [], lb, ank, [], opts);
        % Get the efficiency trajectory:
        effTraj(i,j,:) = EffY(ampV,a(i,j),c(i,j));
    end
end

% Let's obtain the swinging time and switching time:
Ta = (2*ampOptY./a);
Tb = c;
T_ratio = Ta./Tb; % swinging to switching ratio

% % % % Lets obtain the T_beta/T_alpha cost:
% % % cbya = c./(2*ampOptY/a);
% % % 
% % % % Let's now obtain the third surface where efficiency is plotted along the
% % % % other diagonal:
% % % cc = flipud(diag(fliplr(c))); aa = flipud(diag(fliplr(a)));
% % % ccbyaa = flipud(diag(fliplr(cbya))); alphaOpt = flipud(diag(fliplr(ampOptY)));
% % % EOptY = flipud(diag(fliplr(EffOptY)));

% %%
% 
% fSiz = 25;
% 
% % Plot the results:
% f = figure('WindowState','maximized'); % 'WindowState','maximized'
% subplot(3,1,1)
% surf(a,c,ampOptY,'EdgeColor','none');
% hold on;
% ax = gca;
% ax.FontSize = fSiz;
% plot3(aa,cc,alphaOpt,'r','LineWidth',3);
% xlabel('Swing Cost A','fontsize',fSiz); %,'fontweight','bold'
% ylabel('Contact Cost C','fontsize',fSiz); %,'fontweight','bold'
% zlabel('Efficient Swing $\alpha^*$','fontsize',fSiz,'Interpreter','latex'); % ,'fontweight','bold'
% colormap(CUB_sequential); pbaspect([2 2 1]);
% % colorbar();
% 
% subplot(3,1,2)
% surf(a,c,cbya,'EdgeColor','none');
% hold on;
% ax = gca;
% ax.FontSize = fSiz;
% plot3(aa,cc,ccbyaa,'r','LineWidth',3);
% xlabel('A','fontsize',fSiz); %,'fontweight','bold'
% ylabel('C','fontsize',fSiz); %,'fontweight','bold'
% zlabel('C/A','fontsize',fSiz); %,'fontweight','bold'
% colormap(CUB_sequential); pbaspect([2 2 1]);
% % colorbar();
% 
% subplot(3,1,3)
% surf(a,c, 2*sin(alphaOpt)./(c + a.*alphaOpt),'EdgeColor','none');
% hold on;
% ax = gca;
% ax.FontSize = fSiz;
% xlabel('A','fontsize',fSiz); %,'fontweight','bold'
% ylabel('C','fontsize',fSiz); %,'fontweight','bold'
% zlabel('Speed (BL/time)','fontsize',fSiz); %,'fontweight','bold'
% colormap(CUB_sequential); pbaspect([2 2 1]);
% 
% 
% % exportgraphics(f,'OptimalSwing2D_2hybC_nolim.png','Resolution',500);

% %%
% 
% fSiz = 40;
% 
% % Now plot the results along the main diagonal
% figure('WindowState','maximized')
% 
% surf(a,c,ampOptY,'EdgeColor','none');
% hold on;
% ax = gca;
% ax.FontSize = fSiz;
% plot3(aa,cc,alphaOpt,'r','LineWidth',3);
% xlabel('A','fontsize',fSiz); %,'fontweight','bold','Interpreter','latex'
% ylabel('C','fontsize',fSiz); %,'fontweight','bold','Interpreter','latex'
% zlabel('\alpha^*','fontsize',fSiz); % ,'fontweight','bold','Interpreter','latex'
% colormap(CUB_sequential); pbaspect([2 2 1]);
% colorbar();
% 
% % Plot the results using yyaxis
% figure('WindowState','maximized')
% 
% yyaxis left
% plot(cc./aa, 2*sin(alphaOpt)./(cc + aa.*alphaOpt),...
%     'color', circ1, 'LineWidth', 1.2);
% hold on; grid on;
% ax = gca;
% ax.FontSize = fSiz;
% ax.YAxis(1).Color = circ1;
% ylabel('\alpha^*','fontsize',fSiz);
% 
% xlabel('C/A','fontsize',fSiz);
% 
% yyaxis right
% plot(cc./aa, alphaOpt, 'color', circ2, 'LineWidth', 3);
% hold on; grid on;
% ax = gca;
% ax.FontSize = fSiz;
% ax.YAxis(2).Color = circ2;
% ylabel('Speed (BL/time)','fontsize',fSiz);

%%

colorM1 = [255,245,235
254,230,206
253,208,162
253,174,107
253,141,60
241,105,19
217,72,1
166,54,3
127,39,4];

colorM2 = [252,251,253
239,237,245
218,218,235
188,189,220
158,154,200
128,125,186
106,81,163
84,39,143
63,0,125];

fSiz = 40;

% % % OLD STUFF NOT SUPPORTED ANYMORE
% % surf(a,c,ampOptY,'EdgeColor','none');
% % % hold on;
% % ax = gca;
% % ax.FontSize = fSiz;
% % % plot3(aa,cc,alphaOpt,'r','LineWidth',3);
% % xlabel('A','fontsize',fSiz); %,'fontweight','bold','Interpreter','latex'
% % ylabel('C','fontsize',fSiz); %,'fontweight','bold','Interpreter','latex'
% % zlabel('\alpha^*','fontsize',fSiz); % ,'fontweight','bold','Interpreter','latex'
% % colormap(CUB_sequential); pbaspect([3 3 1]);
% % colorbar();

% % az_ang = 35;
% % el_ang = 25;
% % 
% % subplot(1,2,1)
% % surf(c,a,ampOptY,'EdgeColor','none','FaceAlpha',0.4,'FaceColor','k');
% % grid on; hold on;
% % plot3(c(:,1),a(:,1),ampOptY(:,1),'Color','k','LineWidth',3);
% % ax = gca;
% % ax.FontSize = fSiz;
% % % plot3(aa,cc,alphaOpt,'r','LineWidth',3);
% % ylabel('A','fontsize',fSiz); %,'fontweight','bold','Interpreter','latex'
% % xlabel('C','fontsize',fSiz); %,'fontweight','bold','Interpreter','latex'
% % zlabel('\alpha^*','fontsize',fSiz,'Color','k'); % ,'fontweight','bold','Interpreter','latex'
% % ax.ZColor = [0,0,0];
% % pbaspect([3 3 2]); view(az_ang,el_ang);
% % 
% % col2 = [0.6944    0.4340    0.2764]; % copper(100) -- last entry
% % 
% % subplot(1,2,2)
% % surf(c,a,2*sin(ampOptY)./(c + a.*ampOptY),'EdgeColor','none','FaceAlpha'...
% %     ,0.4,'FaceColor',col2);
% % grid on; hold on;
% % plot3(c(:,1),a(:,1),2*sin(ampOptY(:,1))./(c(:,1) + a(:,1).*ampOptY(:,1)),...
% %     'Color',col2,'LineWidth',3);
% % ax = gca;
% % ax.FontSize = fSiz;
% % % plot3(aa,cc,alphaOpt,'r','LineWidth',3);
% % ylabel('A','fontsize',fSiz); %,'fontweight','bold','Interpreter','latex'
% % xlabel('C','fontsize',fSiz); %,'fontweight','bold','Interpreter','latex'
% % zlabel('E_\phi^*','fontsize',fSiz,'Color',col2); % ,'fontweight','bold','Interpreter','latex'
% % ax.ZColor = col2;
% % pbaspect([3 3 2]); view(az_ang,el_ang);
% % 
% % 
% % figure('WindowState','maximized')
% % 
% % pbaspect([3 1 1]);
% % 
% % yyaxis left
% % plot(c(:,1), ampOptY(:,1), 'color', 'k', 'LineWidth', 5);
% % hold on; grid on;
% % ax = gca;
% % ax.FontSize = fSiz+15;
% % ax.YAxis(1).Color = 'k';
% % ylabel('\alpha^*','fontsize',fSiz+15);
% % yticks([1.1,1.2,1.3,1.4,1.5]);
% % 
% % xlim([0.1 1]); xticks([0.2,0.4,0.6,0.8]);
% % 
% % yyaxis right
% % plot(c(:,1), 2*sin(ampOptY(:,1))./(c(:,1) + a(:,1).*ampOptY(:,1)),...
% %     'color', col2, 'LineWidth', 5);
% % hold on; grid on;
% % ax = gca;
% % ax.FontSize = fSiz+15;
% % ax.YAxis(2).Color = col2;
% % ylabel('E_\phi^*','fontsize',fSiz+15);
% % yticks([0,2,4,6,8,10]);
% % 
% % % xlabel('T_\beta','fontsize',fSiz+15);

% Plot the results
figure('WindowState','maximized')

az_ang = 35;
el_ang = 25;

subplot(1,2,1) % lets use the last slice instead for our comparison
surf(Tb,Ta,ampOptY,'EdgeColor','none','FaceAlpha',0.4,'FaceColor','k');
grid on; hold on;
plot3(Tb(:,end),Ta(:,end),ampOptY(:,end),'Color','k','LineWidth',3);
ax = gca;
ax.FontSize = fSiz;
% plot3(aa,cc,alphaOpt,'r','LineWidth',3);
ylabel('$T_\alpha$','fontsize',fSiz,'Interpreter','latex'); %,'fontweight','bold','Interpreter','latex'
xlabel('$T_\beta$','fontsize',fSiz,'Interpreter','latex'); %,'fontweight','bold','Interpreter','latex'
zlabel('$\hat{\alpha}^*$','fontsize',fSiz,'Interpreter','latex'); % ,'fontweight','bold','Interpreter','latex'
ax.ZColor = [0,0,0];
pbaspect([3 3 2]); view(az_ang,el_ang);

col2 = [0.6944    0.4340    0.2764]; % copper(100) -- last entry

subplot(1,2,2) % negative sign added since it is solved as a minimization problem
surf(Tb,Ta,-EffOptY,'EdgeColor','none','FaceAlpha',0.4,'FaceColor',col2);
grid on; hold on;
plot3(Tb(:,end),Ta(:,end),-EffOptY(:,end),'Color',col2,'LineWidth',3);
ax = gca;
ax.FontSize = fSiz;
% plot3(aa,cc,alphaOpt,'r','LineWidth',3);
ylabel('$T_\alpha$','fontsize',fSiz,'Interpreter','latex'); %,'fontweight','bold','Interpreter','latex'
xlabel('$T_\beta$','fontsize',fSiz,'Interpreter','latex'); %,'fontweight','bold','Interpreter','latex'
zlabel('$E_{\phi}^*$','fontsize',fSiz,'Color',col2,'Interpreter','latex'); % ,'fontweight','bold','Interpreter','latex'
ax.ZColor = col2;
pbaspect([3 3 2]); view(az_ang,el_ang);


figure('WindowState','maximized')

pbaspect([3 1 1]);

yyaxis left
plot(T_ratio(:,end), ampOptY(:,end), 'color', 'k', 'LineWidth', 5);
hold on; grid on;
ax = gca;
ax.FontSize = fSiz+15;
ax.YAxis(1).Color = 'k';
ax.TickLabelInterpreter='latex';
ylabel('$\hat{\alpha}^*$','fontsize',fSiz+15,'Interpreter','latex','Rotation',0);
yticks([0.8,1,1.2]);

xlim([min(T_ratio(:,end)) max(T_ratio(:,end))]); 
xticks([1,2,3,4]);
xlabel('$\frac{T_{\alpha}}{T_{\beta}}$','fontsize',fSiz+15,...
    'Interpreter','latex');

yyaxis right
plot(T_ratio(:,end), -EffOptY(:,end),...
    'color', col2, 'LineWidth', 5);
hold on; grid on;
ax = gca;
ax.FontSize = fSiz+15;
ax.YAxis(2).Color = col2;
ax.TickLabelInterpreter='latex';
ylabel('$E_{\phi}^y$','fontsize',fSiz+15,'Interpreter','latex','Rotation',0,'HorizontalAlignment','left');
yticks([1.2,1.8,2.4]);

% axis tight;

% xlabel('T_\beta','fontsize',fSiz+15);

% % % % %% OLD STUFF NOT SUPPORTED ANYMORE
% % % % figure('WindowState','maximized'); % f = 
% % % % surf(a,c,cbya,'EdgeColor','none');
% % % % hold on;
% % % % ax = gca;
% % % % ax.FontSize = fSiz;
% % % % plot3(aa,cc,ccbyaa,'r','LineWidth',3);
% % % % xlabel('A','fontweight','bold','fontsize',fSiz);
% % % % ylabel('C','fontweight','bold','fontsize',fSiz);
% % % % % zlabel('C/A','fontweight','bold','fontsize',fSiz);
% % % % text(1.015,1.03,'C/A','fontweight','bold','fontsize',fSiz);
% % % % colormap(CUB_sequential); pbaspect([1 1 1]); colorbar;
% % % % view(0,90);
% % % % % % % % % 
% % % % % % % % % exportgraphics(f,'OptimalSwing2D_2hybC_nolim_justCbyA_Overleaf.png','Resolution',500);

%%

% % idx = logical(repmat(fliplr(eye(size(effTraj,1))),1,1,size(effTraj,3)));
% % EffY_sweep = -effTraj(idx); EffY_sweep = reshape(EffY_sweep,size(a));
% % 
% % ca = nan(size(ccbyaa)); e = nan(size(ca)); aaa = nan(size(e));
% % for i = 1:length(ampV) % obtain the max efficiency and corresponding c/a values:
% %     [eM, idx_eM] = max(EffY_sweep(i,:));
% %     ca(i) = ccbyaa(idx_eM);
% %     aaa(i) = aa(idx_eM);
% %     e(i) = eM;
% % end
% % 
% % % Since the c/a variation is a little too squiggly, let's use a moving
% % % mean.
% % ca = smooth(ca);
% % 
% % mulC = 1/255;
% % temp1 = mulC*[178,171,210]; temp2 = mulC*([128,115,172]); temp3 = mulC*[84,39,136];
% % 
% % fSiz = 25;
% % 
% % f = figure('WindowState','Maximize');
% % 
% % subplot(2,1,1);
% % surf(log10(ccbyaa),ampV,EffY_sweep,'EdgeColor','none','FaceAlpha',1.0);
% % ax = gca;
% % ax.FontSize = fSiz;
% % hold on;
% % plot3(log10(ca),ampV,e,'Color','k','LineWidth',3);
% % 
% % plot3(log10(ccbyaa(1))*ones(size(ampV)),ampV,EffY_sweep(:,1),'Color',temp3,'LineWidth',5);
% % plot3(log10(ccbyaa(7))*ones(size(ampV)),ampV,EffY_sweep(:,7),'Color',temp2,'LineWidth',4);
% % % plot3(log10(ccbyaa(50))*ones(size(ampV)),ampV,EffY_sweep(:,50),'Color',temp2,'LineWidth',4);
% % plot3(log10(ccbyaa(end))*ones(size(ampV)),ampV,EffY_sweep(:,end),'Color',temp1,'LineWidth',3);
% % 
% % 
% % 
% % text(log10(0.8),1,0.9,'\epsilon_\phi^y','fontweight','bold','fontsize',fSiz,'Color','k');
% % 
% % xlabel('log_{10} (C/A)','fontweight','bold','fontsize',fSiz);
% % ylabel('\alpha_0 (rad)','fontweight','bold','fontsize',fSiz);
% % zlabel('E_\phi^y','fontweight','bold','fontsize',fSiz);
% % yticks([0, 1, pi/2]); yticklabels({'0', '1', '\pi/2'}); %\pi/2
% % xticks(linspace(1,10,10));
% % % view(-20,25);
% % view(-45,25);
% % colormap(CUB); 
% % colorbar();
% % pbaspect([2 1 0.5]);
% % 
% % s2 = subplot(2,1,2);
% % surf(log10(ccbyaa),ampV,EffY_sweep,'EdgeColor','none','FaceAlpha',0.2);
% % % set(s2, 'Color', mulC*[217,217,217]);
% % ax = gca;
% % ax.FontSize = fSiz;
% % hold on;
% % 
% % plot3(ca,ampV,e,'Color','k','LineWidth',4);
% % 
% % 
% % 
% % grid off;
% % 
% % plot3(log10(ccbyaa(1))*ones(size(ampV)),ampV,EffY_sweep(:,1),'Color',temp3,'LineWidth',5);
% % plot3(log10(ccbyaa(7))*ones(size(ampV)),ampV,EffY_sweep(:,7),'Color',temp2,'LineWidth',4);
% % % plot3(log10(ccbyaa(50))*ones(size(ampV)),ampV,EffY_sweep(:,50),'Color',temp2,'LineWidth',4);
% % plot3(log10(ccbyaa(end))*ones(size(ampV)),ampV,EffY_sweep(:,end),'Color',temp1,'LineWidth',3);
% % 
% % text(-1,1.02,0.85,'\epsilon_\phi^y','fontweight','bold','fontsize',fSiz-5,'Color','k');
% % 
% % text(0.8,1,0.08,num2str(round(ccbyaa(1),1)),'fontweight','bold','fontsize',fSiz-5,'Color',temp3);
% % text(0.8,1,0.26,num2str(round(ccbyaa(7),1)),'fontweight','bold','fontsize',fSiz-5,'Color',temp2);
% % % text(0.8,1,0.26,num2str(round(ccbyaa(50),1)),'fontweight','bold','fontsize',fSiz-5,'Color',temp2);
% % text(0.8,1,0.63,['C/A = ',num2str(round(ccbyaa(end),1))],'fontweight','bold','fontsize',fSiz-5,'Color',temp1);
% % 
% % % xlabel('C/A','fontweight','bold','fontsize',fSiz);
% % ylabel('\alpha_0 (rad)','fontweight','bold','fontsize',fSiz);
% % zlabel('E_\phi^y','fontweight','bold','fontsize',fSiz);
% % yticks([0, 1, pi/2]); yticklabels({'0', '1', '\pi/2'}); %\pi/2
% % % xticks(linspace(1,10,10));
% % view(-90,0);
% % colormap(CUB); 
% % % colorbar();
% % pbaspect([2 1 0.5]);

%%

% CUB = flipud([158,202,225;
%             107,174,214;
%             66,146,198;
%             33,113,181;
%             8,81,156;
%             8,48,107])*(1/255);
% 
% CUB = interp1(linspace(0,discN,size(CUB,1)), CUB, linspace(0,discN,251));
% 
% % Set the font size:
% fSiz = 25;
% 
% % %%
% % A random optimal efficiency value
% r1idx = randi(length(A)); r2idx = randi(length(C));
% 
% % Plot the efficiency as a function of swing
% f = figure('WindowState','maximized');
% 
% subplot(1,2,1);
% plot(ampV,reshape(-effTraj(r1idx,r2idx,:),1,size(effTraj,3)), 'k', 'linewidth', 3);
% hold on;
% 
% ax = gca;
% ax.FontSize = fSiz; 
% 
% xline(ampOptY(r1idx,r2idx), 'k', 'linewidth', 1.5); text(1,0.2,'Optimal Swing','FontSize',fSiz,'Color','k');
% yline(-EffY(ampOptY(r1idx,r2idx),a(r1idx,r2idx),c(r1idx,r2idx)), 'k', 'linewidth', 1.5);
% text(2,1.34,'Peak Efficiency','FontSize',fSiz,'Color','k');
% 
% xlabel('Swing Amplitude \alpha_0 range (rad)','FontSize',fSiz);
% ylabel('Efficiency E_\phi','FontSize',fSiz); % defined in a negative sense and inverted
% 
% ylim([0 1.5]);
% 
% subplot(1,2,2);
% 
% yyaxis left
% plot(ampV,c(r1idx,r2idx)./(a(r1idx,r2idx)*ampV), 'b', 'linewidth', 3);
% hold on;
% 
% xline(ampOptY(r1idx,r2idx), 'k', 'linewidth', 1.5);
% yline(c(r1idx,r2idx)./(a(r1idx,r2idx)*ampOptY(r1idx,r2idx)), 'k', 'linewidth', 1.5);
% text(1.5,0.4,'Optimal Ratio','FontSize',fSiz,'Color','k');
% % yline(-EffY(ampOptY(r1idx,r2idx),a(r1idx,r2idx),c(r1idx,r2idx)), 'r', 'linewidth', 1.5);
% 
% text(2,2.75,['Contact Cost = ', num2str(round(c(r1idx,r2idx),2))],'FontSize',fSiz);
% text(2,2.5,['Swing Cost = ', num2str(round(a(r1idx,r2idx),2))],'FontSize',fSiz);
% 
% ax = gca;
% ax.FontSize = fSiz;
% 
% % xlabel('Swing Amplitude Range \alpha_0 (rad)','FontSize',fSiz);
% ylabel('Shape-change Cost Ratio ( C / A\alpha_0 )','FontSize',fSiz); % defined in a negative sense and inverted
% axis([-0.25 3.5 0 3]);
% 
% yyaxis right
% plot(ampV,sin(ampV)./(a(r1idx,r2idx)*ampV), 'r--', 'linewidth', 3);
% hold on;
% ax = gca;
% ax.FontSize = fSiz;
% ylabel('E^\infty_\phi = sin\alpha_0 / \alpha_0 ','FontSize',fSiz);
% yline(sin(ampOptY(r1idx,r2idx))/(a(r1idx,r2idx)*ampOptY(r1idx,r2idx)), 'k', 'linewidth', 1.5);

% exportgraphics(f,'OptEffTrajRandPt_2hybC_nolim.png','Resolution',500);

% %% OPTIMAL Y-DISPLACEMENT SWING AMPLITUDE FOR 2-CONTACT SYSTEM
% 
% % clear the variables to prevent overwriting during multiple runs:
% clear gcirc_small gcirc_large
% 
% % Let's take the variation along the other diagonal to highlight our point
% % on contact bandwidth to swing bandwidth,
% effS = fliplr(diag(fliplr(ampOptS))); %effY = effY(1:54);
% xSwing = fliplr(A);
% yContact = C;
% cost = yContact./xSwing; %cost = cost(1:54);
% 
% % Let's pick 2 points to show the trajectories for below
% cost_2plot = [cost(1) cost(end)];
% swing_amplitude_2plot = [effS(1) effS(end)];
% % SET IT UP IN THE GUI in such a way that we evaluate each gait
% % sequentially. Look below for further instructions.
% 
% % Open the GUI load the system and shape change, plot the trajectory and
% % BVI first. Then go into the
% % GeometricSystemPlotter\UserFiles\GenericUser\sysplotter_data to load the
% % appropriate files. You need to do this twice for each gait we are trying
% % to plot.
% 
% % % Lets load the calculation files from sysplotter - p structure for the
% % % larger displacement stride
% % load('sysf_hybridC_2c1r__shchf_hybridC_Square'); gcirc_large = p.G_locus_full{1,1}.G(:,1:2);
% % 
% % %%
% % % p structure for the smaller displacement stride -- MAKE CHANGES TO THE
% % % SHAPE CHANGE FILE FIRST AND RUN THE NEW PLOT TO CHECK
% % load('sysf_hybridC_2c1r__shchf_hybridC_Square'); gcirc_small = p.G_locus_full{1,1}.G(:,1:2);
% 
% % % Let's save the data we have locally:
% % save('OptimalYDisplacementSwing_2hybC.mat',"gcirc_small","gcirc_large");
% 
% % Let's plot this out:
% f = figure(); % 'WindowState','maximized'
% 
% % Let's make two colors for the small and large gaits:
% cR = [215,48,39]/255; cB = [69,117,180]/255;
% 
% subplot(2,1,1)
% plot(effS, cost, 'k', 'linewidth', 3,'DisplayName','Optimal Swing Amplitude');
% grid on; hold on;
% plot(effS(1), cost(1), 'color',cB,'Marker','square', 'markersize', 25,'linewidth',2,...
%     'DisplayName','Optimal Swing Amplitude');
% plot(effS(end), cost(end), 'color',cR,'Marker','square', 'markersize', 25,'linewidth',4,...
%     'DisplayName','Optimal Swing Amplitude');
% % yline(swing_amplitude_2plot(1),'color',[178, 24, 43]/255,'LineStyle','--','linewidth',2);
% % yline(swing_amplitude_2plot(2),'color',[178, 24, 43]/255,'LineStyle','--','linewidth',4);
% % xline(cost_2plot(1),'color',[178, 24, 43]/255,'LineStyle','--','linewidth',2);
% % xline(cost_2plot(2),'color',[178, 24, 43]/255,'LineStyle','--','linewidth',4);
% ylabel('Cost Ratio (C/A)','fontweight','bold','fontsize',fontS);
% xlabel('Angular swing amplitude (rad)','fontweight','bold','fontsize',fontS);
% title('Efficient swing for rate-limited 2-contact system');
% ylim([-0.1 11]); 
% % axis equal;
% 
% % Load the data we have for the optimal gait displacements:
% load('OptimalYDisplacementSwing_2hybC.mat');
% 
% % Offset the large displacement to match the optimal swing amplitude plot
% % xlines to show correspondence:
% off_small = -0.5; off_large = 0.5;
% gcirc_large(:,1) = gcirc_large(:,1) + off_large;
% gcirc_small(:,1) = gcirc_small(:,1) + off_small;
% 
% % Robot body-length (BL) is defined as twice the leg-length for the 2
% % contact case:
% bl = 2/3; % this choice helped fit the animation properly, nothing more.
% gcirc_small = gcirc_small/bl; gcirc_large = gcirc_large/bl;
% 
% % Make the second subplot:
% subplot(2,1,2)
% plot(gcirc_small(:,1),gcirc_small(:,2),'color',cB,'linewidth',2);
% grid on; hold on;
% yline(0,'k--','LineWidth',1.2); 
% xline(off_small/bl,'k--','LineWidth',1.2); xline(off_large/bl,'k--','LineWidth',1.2);
% plot(gcirc_large(:,1),gcirc_large(:,2),'color',cR,'linewidth',4);
% axis equal; ylim([-0.8 0.1]); xticks([]); 
% ylabel('Displacement (BLs)','fontweight','bold','fontsize',fontS);
% 
% % MAKE SURE YOU'RE IN YOUR WORK FOLDER NOW: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % exportgraphics(f,'OptimalSwingTrajectoryPlusVar_2hybC.png','Resolution',500);
% 
% %%% Method to switch off annotation
% % h = plot(x, sin(x));
% % h.Annotation.LegendInformation.IconDisplayStyle = 'off';

%% A PANELS APPROACH: EXPRESSING THE VERTICAL DISTANCE BETWEEN HYBRID CONTACT VECTOR FIELD
% We make use of analytical results from the 2-contact and the 3-contact
% problems. We can think of these veritcal "panels" as average CCFs.

fSiz = 40;
discN = 100;
% CUB = flipud([215,48,39;
%         244,109,67;
%         253,174,97;
%         254,224,144;
%         255,255,191;
%         224,243,248;
%         171,217,233;
%         116,173,209;
%         69,117,180])*(1/255);
% 
% CUB = interp1(linspace(0,100,size(CUB,1)), CUB, linspace(0,100,251));

% Function that computes the connection
Ax = @(a) sin(a); Ay = @(a) -cos(a);

% Compute the panels for the two contact system:
alpha_dom = [-pi/2-pi/36,pi/2+pi/36]; % leg swing domain - continuous 95degs
contact_dom = [-1.1,1.1]; % contact domain - discrete (slight padding added)

% Create a grid:
[a,c] = meshgrid(linspace(alpha_dom(1),alpha_dom(2),discN),...
    linspace(contact_dom(1),contact_dom(2),discN));

% Compute the vertical panel value:
panel_2C_X = Ax(a) - Ax(a + pi);
panel_2C_Y = Ay(a) - Ay(a + pi);

%%
% Plot the panels:
f = figure('WindowState','maximized'); %'WindowState','maximized'

fSiz = 70;

subplot(1,2,1)
surf(a,c,panel_2C_X,'EdgeColor','none');
ax = gca;
pbaspect([2,1,1]);
ax.FontSize = fSiz;
view(0,90);
% ylabel('contact state','fontweight','bold','fontsize',fSiz,...
%     'VerticalAlignment','middle'); %'Rotation',0,
% ylabel('$c$','fontsize',fSiz,...
%     'VerticalAlignment','middle','Interpreter','latex'); %'Rotation',0,'fontweight','bold',
yticks([-1 1]);
% yticks([]);
% yticklabels({'\{1,0\}','\{0,1\}'});
% xlabel('\alpha','fontweight','bold','fontsize',fSiz);
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
% xticks([]);
xtickangle(0);
% title('\Delta X','fontweight','bold','fontsize',fSiz); % axis square;
caxis([-2 2]);
colormap(CUB); colorbar('Ticks',linspace(-2,2,3),'FontSize',fSiz);
% colorbar('Ticks',[]);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.1 1.1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)
surf(a,c,panel_2C_Y,'EdgeColor','none');
ax = gca;
pbaspect([2,1,1]);
ax.FontSize = fSiz;
view(0,90);
% ylabel('$c$','fontsize',fSiz,...
%     'VerticalAlignment','middle','Interpreter','latex'); %'Rotation',0,'fontweight','bold',
yticks([-1 1]);
% yticks([]);
% yticklabels({'\{1,0\}','\{0,1\}'});
% xlabel('\alpha','fontweight','bold','fontsize',fSiz);
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0);
% xticks([]);
% title('\Delta Y','fontweight','bold','fontsize',fSiz); % axis square;
caxis([-2 2]);
colormap(CUB); 
colorbar('Ticks',linspace(-2,2,3),'FontSize',fSiz);
% colorbar('Ticks',[]);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.1 1.1]);

%%

fSiz = 40;

figure() %'WindowState','maximized'
surf(a,c,panel_2C_Y,'EdgeColor','none');
ax = gca;
pbaspect([2,1,1]);
ax.FontSize = fSiz;
view(0,90);
ylabel('$c$','fontsize',fSiz,...
    'VerticalAlignment','middle','Interpreter','latex'); %'Rotation',0,'fontweight','bold',
yticks([-1 1]);
% yticklabels({'\{1,0\}','\{0,1\}'});
xlabel('\alpha','fontweight','bold','fontsize',fSiz);
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
% title('\Delta Y','fontweight','bold','fontsize',fSiz); % axis square;
caxis([-2 2]);
colormap(CUB); colorbar('Ticks',linspace(-2,2,5),'FontSize',fSiz); 
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.1 1.1]);

%% NEW PANELS 2C - 1D
f = figure('WindowState','maximized');
subplot(1,2,1)
surf(a,c,panel_2C_X,'EdgeColor','none');
ax = gca;
pbaspect([10,1,1]);
ax.FontSize = fSiz;
view(0,90);
yticks([]);
xticks([-pi/2 pi/2]); xticklabels({'-\pi/2','\pi/2'});
xtickangle(0);
caxis([-2 2]);
colormap(CUB); 
colorbar('Ticks',linspace(-2,2,2),'FontSize',fSiz); 
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.1 1.1]);

subplot(1,2,2)
surf(a,c,panel_2C_Y,'EdgeColor','none');
ax = gca;
pbaspect([10,1,1]);
ax.FontSize = fSiz;
view(0,90);
yticks([]);
xticks([-pi/2 pi/2]); xticklabels({'-\pi/2','\pi/2'});
xtickangle(0);
caxis([-2 2]);
colormap(CUB); 
colorbar('Ticks',linspace(-2,2,2),'FontSize',fSiz); 
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.1 1.1]);



%%

% % Plot the panels - Vertical
% f = figure('WindowState','maximized'); %'WindowState','maximized'
% 
% subplot(2,1,1)
% surf(a,c,panel_2C_X,'EdgeColor','none');
% ax = gca;
% ax.FontSize = fSiz;
% view(0,90);
% ylabel('contact state','fontweight','bold','fontsize',fSiz,...
%     'VerticalAlignment','middle'); %'Rotation',0,
% yticks([-1 1]);
% yticklabels({'\{1,0\}','\{0,1\}'});
% xlabel('\alpha (rad)','fontweight','bold','fontsize',fSiz);
% xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','pi/2'});
% axis square; title('\Delta X','fontweight','bold','fontsize',fSiz);
% colormap(CUB); colorbar; xlim([-pi/2 pi/2]);
% 
% subplot(2,1,2)
% surf(a,c,panel_2C_Y,'EdgeColor','none');
% ax = gca;
% ax.FontSize = fSiz;
% view(0,90);
% yticks([-1 1]);
% yticklabels({'\{1,0\}','\{0,1\}'});
% xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','pi/2'});
% axis square; title('\Delta Y','fontweight','bold','fontsize',fSiz);
% colormap(CUB); colorbar; xlim([-pi/2 pi/2]);

% sgtitle('Hybrid 2-Contact System Avg CCF','fontweight','bold','fontsize',fontS+2);

% MAKE SURE YOU'RE IN YOUR WORK FOLDER NOW: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% exportgraphics(f,'panels_2C.png','Resolution',500);

% % % % %%
% % % % 
% % % % fSiz = 20;
% % % % 
% % % % % Compute the panesl for each contact switch (3 switches):
% % % % panel12_3C_X = Ax(a) - Ax(a + 2*pi/3);
% % % % panel12_3C_Y = Ay(a) - Ay(a + 2*pi/3);
% % % % 
% % % % panel23_3C_X = Ax(a + 2*pi/3) - Ax(a + 4*pi/3);
% % % % panel23_3C_Y = Ay(a + 2*pi/3) - Ay(a + 4*pi/3);
% % % % 
% % % % panel31_3C_X = Ax(a + 4*pi/3) - Ax(a);
% % % % panel31_3C_Y = Ay(a + 4*pi/3) - Ay(a);
% % % % 
% % % % % Plot the panels in a CYCLIC order:
% % % % f = figure('WindowState','maximized'); %'WindowState','maximized'
% % % % 
% % % % subplot(3,2,1)
% % % % surf(a,c,panel12_3C_X,'EdgeColor','none');
% % % % view(0,90);
% % % % ax = gca;
% % % % ax.FontSize = fSiz;
% % % % % ylabel('c','fontweight','bold','fontsize',fSiz,...
% % % % %     'Rotation',0,'VerticalAlignment','middle');
% % % % yticks([-1 1]);
% % % % yticklabels({'\{1,0,0\}','\{0,1,0\}'});
% % % % xlabel('\alpha (rad)','fontweight','bold','fontsize',fSiz);
% % % % xticks([-1 -0.5 0 0.5 1]); 
% % % % axis square; title('\Delta X','fontweight','bold','fontsize',fSiz);
% % % % colormap(CUB); colorbar; xlim([-1 1]);
% % % % 
% % % % subplot(3,2,2)
% % % % surf(a,c,panel12_3C_Y,'EdgeColor','none');
% % % % view(0,90);
% % % % ax = gca;
% % % % ax.FontSize = fSiz;
% % % % yticks([-1 1]);
% % % % yticklabels({'\{1,0,0\}','\{0,1,0\}'});
% % % % xticks([-1 -0.5 0 0.5 1]); 
% % % % axis square; title('\Delta Y','fontweight','bold','fontsize',fSiz);
% % % % colormap(CUB); colorbar; xlim([-1 1]);
% % % % 
% % % % subplot(3,2,3)
% % % % surf(a,c,panel23_3C_X,'EdgeColor','none');
% % % % view(0,90);
% % % % ax = gca;
% % % % ax.FontSize = fSiz;
% % % % yticks([-1 1]);
% % % % yticklabels({'\{0,1,0\}','\{0,0,1\}'});
% % % % xticks([-1 -0.5 0 0.5 1]); 
% % % % axis square; colormap(CUB); colorbar; xlim([-1 1]);
% % % % 
% % % % subplot(3,2,4)
% % % % surf(a,c,panel23_3C_Y,'EdgeColor','none');
% % % % view(0,90);
% % % % ax = gca;
% % % % ax.FontSize = fSiz;
% % % % yticks([-1 1]);
% % % % yticklabels({'\{0,1,0\}','\{0,0,1\}'});
% % % % xticks([-1 -0.5 0 0.5 1]); 
% % % % axis square; colormap(CUB); colorbar; xlim([-1 1]);
% % % % 
% % % % subplot(3,2,5)
% % % % surf(a,c,panel31_3C_X,'EdgeColor','none');
% % % % view(0,90);
% % % % ax = gca;
% % % % ax.FontSize = fSiz;
% % % % yticks([-1 1]);
% % % % yticklabels({'\{0,0,1\}','\{1,0,0\}'});
% % % % xticks([-1 -0.5 0 0.5 1]); 
% % % % axis square; colormap(CUB); colorbar; xlim([-1 1]);
% % % % 
% % % % subplot(3,2,6)
% % % % surf(a,c,panel31_3C_Y,'EdgeColor','none');
% % % % view(0,90);
% % % % ax = gca;
% % % % ax.FontSize = fSiz;
% % % % yticks([-1 1]);
% % % % yticklabels({'\{0,0,1\}','\{1,0,0\}'});
% % % % xticks([-1 -0.5 0 0.5 1]); 
% % % % axis square; colormap(CUB); colorbar; xlim([-1 1]);
% % % % 
% % % % % sgtitle('Hybrid 3-Contact System Avg CCF','fontweight','bold','fontsize',fSiz+2);
% % % % 
% % % % % MAKE SURE YOU'RE IN YOUR WORK FOLDER NOW: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % % % % exportgraphics(f,'panels_3C.png','Resolution',500);
% % % % 
% % % % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%

fSiz = 50; discN = 100;

[~,c13] = meshgrid(linspace(alpha_dom(1),alpha_dom(2),discN),...
    linspace(contact_dom(1),contact_dom(2),discN));
[~,c12] = meshgrid(linspace(alpha_dom(1),alpha_dom(2),discN),...
    linspace(contact_dom(1),0.05,discN));
[~,c23] = meshgrid(linspace(alpha_dom(1),alpha_dom(2),discN),...
    linspace(-0.05,contact_dom(2),discN));

% Compute the panesl for each contact switch (3 switches):
panel12_3C_X = Ax(a) - Ax(a + 2*pi/3);
panel12_3C_Y = Ay(a) - Ay(a + 2*pi/3);

panel23_3C_X = Ax(a + 2*pi/3) - Ax(a + 4*pi/3);
panel23_3C_Y = Ay(a + 2*pi/3) - Ay(a + 4*pi/3);

panel31_3C_X = -(Ax(a + 4*pi/3) - Ax(a)); % added a minus sign for the paper now panel 13
panel31_3C_Y = -(Ay(a + 4*pi/3) - Ay(a));

% Plot the panels in a CYCLIC order:
f = figure('WindowState','maximized'); %'WindowState','maximized'

val = 1.28;

if max(panel23_3C_X,[],'all') > eps
    maxZx = floor(10*max(panel23_3C_X,[],'all'))*0.1;
else
    maxZx = ceil(10*max(panel23_3C_X,[],'all'))*0.1;
end
if min(panel23_3C_X,[],'all') < -eps
    minZx = ceil(10*min(panel23_3C_X,[],'all'))*0.1;
else
    minZx = floor(10*min(panel23_3C_X,[],'all'))*0.1;
end

subplot(3,2,3)
surf(a,c23,panel23_3C_X,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
% y2 = ylabel('p_{23}','fontweight','bold','fontsize',fSiz,'Rotation',0,...
%     'HorizontalAlignment','center');
% y2.Position(1) = y2.Position(1) - val;
% ylabel('c','fontweight','bold','fontsize',fSiz,...
%     'Rotation',0,'VerticalAlignment','middle');
yticks([0 1]);
% yticklabels({'\{1,0,0\}','\{0,1,0\}'});
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0);
% axis square;
pbaspect([3,1,1]);
colormap(CUB); 
% colorbar('Ticks',linspace(minZx,maxZx,3),'FontSize',fSiz);
colorbar('Ticks',linspace(-2,2,5),'FontSize',fSiz);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-0.05 1.05]);

if max(panel23_3C_Y,[],'all') > eps
    maxZy = floor(10*max(panel23_3C_Y,[],'all'))*0.1;
else
    maxZy = ceil(10*max(panel23_3C_Y,[],'all'))*0.1;
end
if min(panel23_3C_Y,[],'all') < -eps
    minZy = ceil(10*min(panel23_3C_Y,[],'all'))*0.1;
else
    minZy = floor(10*min(panel23_3C_Y,[],'all'))*0.1;
end

subplot(3,2,4)
surf(a,c23,panel23_3C_Y,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
yticks([0 1]);
% yticklabels({'\{1,0,0\}','\{0,1,0\}'});
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0);
% axis square; 
pbaspect([3,1,1]);
colormap(CUB); 
% colorbar('Ticks',linspace(minZy,maxZy,3),'FontSize',fSiz);
colorbar('Ticks',linspace(-2,2,5),'FontSize',fSiz);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-0.05 1.05]);

if max(panel12_3C_X,[],'all') > eps
    maxZx = floor(10*max(panel12_3C_X,[],'all'))*0.1;
else
    maxZx = ceil(10*max(panel12_3C_X,[],'all'))*0.1;
end
if min(panel12_3C_X,[],'all') < -eps
    minZx = ceil(10*min(panel12_3C_X,[],'all'))*0.1;
else
    minZx = floor(10*min(panel12_3C_X,[],'all'))*0.1;
end

subplot(3,2,1)
surf(a,c12,panel12_3C_X,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
y1 = ylabel('$c$','fontsize',fSiz,...
    'HorizontalAlignment','center','Interpreter','latex'); % ,'fontweight','bold','Rotation',0,
y1.Position(1) = y1.Position(1) + val;
xlabel('\alpha','fontweight','bold','fontsize',fSiz);
yticks([-1 0]);
% yticklabels({'\{0,1,0\}','\{0,0,1\}'});
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0);
% title('\Delta X','fontweight','bold','fontsize',fSiz);
% axis square; 
pbaspect([3,1,1]);
colormap(CUB); 
% colorbar('Ticks',linspace(minZx,maxZx,3),'FontSize',fSiz);
colorbar('Ticks',linspace(-2,2,5),'FontSize',fSiz);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.05 0.05]);

if max(panel12_3C_Y,[],'all') > eps
    maxZy = floor(10*max(panel12_3C_Y,[],'all'))*0.1;
else
    maxZy = ceil(10*max(panel12_3C_Y,[],'all'))*0.1;
end
if min(panel12_3C_Y,[],'all') < -eps
    minZy = ceil(10*min(panel12_3C_Y,[],'all'))*0.1;
else
    minZy = floor(10*min(panel12_3C_Y,[],'all'))*0.1;
end

subplot(3,2,2)
surf(a,c12,panel12_3C_Y,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
yticks([-1 0]);
% yticklabels({'\{0,1,0\}','\{0,0,1\}'});
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0);
% title('\Delta Y','fontweight','bold','fontsize',fSiz);
% axis square; 
pbaspect([3,1,1]);
colormap(CUB); 
% colorbar('Ticks',linspace(minZy,maxZy,3),'FontSize',fSiz);
colorbar('Ticks',linspace(-2,2,5),'FontSize',fSiz);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.05 0.05]);

if max(panel31_3C_X,[],'all') > eps
    maxZx = floor(10*max(panel31_3C_X,[],'all'))*0.1;
else
    maxZx = ceil(10*max(panel31_3C_X,[],'all'))*0.1;
end
if min(panel31_3C_X,[],'all') < -eps
    minZx = ceil(10*min(panel31_3C_X,[],'all'))*0.1;
else
    minZx = floor(10*min(panel31_3C_X,[],'all'))*0.1;
end

subplot(3,2,5)
surf(a,c13,panel31_3C_X,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
% y3 = ylabel('p_{23}','fontweight','bold','fontsize',fSiz,'Rotation',0,...
%     'HorizontalAlignment','center');
% y3.Position(1) = y3.Position(1) - val;
yticks([-1 1]);
% yticklabels({'\{1,0,0\}','\{0,0,1\}'});
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0);
% axis square; 
pbaspect([3,1,1]);
colormap(CUB); 
% colorbar('Ticks',linspace(minZx,maxZx,3),'FontSize',fSiz);
colorbar('Ticks',linspace(-2,2,5),'FontSize',fSiz);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.1 1.1]);

if max(panel31_3C_Y,[],'all') > eps
    maxZy = floor(10*max(panel31_3C_Y,[],'all'))*0.1;
else
    maxZy = ceil(10*max(panel31_3C_Y,[],'all'))*0.1;
end
if min(panel31_3C_Y,[],'all') < -eps
    minZy = ceil(10*min(panel31_3C_Y,[],'all'))*0.1;
else
    minZy = floor(10*min(panel31_3C_Y,[],'all'))*0.1;
end

subplot(3,2,6)
surf(a,c13,panel31_3C_Y,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
yticks([-1 1]);
% yticklabels({'\{1,0,0\}','\{0,0,1\}'});
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0);
% axis square; 
pbaspect([3,1,1]);
colormap(CUB); 
% colorbar('Ticks',linspace(minZy,maxZy,3),'FontSize',fSiz);
colorbar('Ticks',linspace(-2,2,5),'FontSize',fSiz);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.1 1.1]);

% sgtitle('Hybrid 3-Contact System Avg CCF','fontweight','bold','fontsize',fSiz+2);

% MAKE SURE YOU'RE IN YOUR WORK FOLDER NOW: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% exportgraphics(f,'panels_3C_matched.png','Resolution',500);

%% NEW PANELS 3C

fSiz = 50; discN = 100;

[~,c13] = meshgrid(linspace(alpha_dom(1),alpha_dom(2),discN),...
    linspace(contact_dom(1),contact_dom(2),discN));
[~,c12] = meshgrid(linspace(alpha_dom(1),alpha_dom(2),discN),...
    linspace(contact_dom(1),0.05,discN));
[~,c23] = meshgrid(linspace(alpha_dom(1),alpha_dom(2),discN),...
    linspace(-0.05,contact_dom(2),discN));

panel12_3C_X = Ax(a) - Ax(a + 2*pi/3);
panel12_3C_Y = Ay(a) - Ay(a + 2*pi/3);

panel23_3C_X = Ax(a + 2*pi/3) - Ax(a + 4*pi/3);
panel23_3C_Y = Ay(a + 2*pi/3) - Ay(a + 4*pi/3);

panel31_3C_X = -(Ax(a + 4*pi/3) - Ax(a));
panel31_3C_Y = -(Ay(a + 4*pi/3) - Ay(a));

f = figure('WindowState','maximized');

val = 1.28;

if max(panel23_3C_X,[],'all') > eps
    maxZx = floor(10*max(panel23_3C_X,[],'all'))*0.1;
else
    maxZx = ceil(10*max(panel23_3C_X,[],'all'))*0.1;
end
if min(panel23_3C_X,[],'all') < -eps
    minZx = ceil(10*min(panel23_3C_X,[],'all'))*0.1;
else
    minZx = floor(10*min(panel23_3C_X,[],'all'))*0.1;
end

subplot(3,2,3)
surf(a,c23,panel23_3C_X,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
yticks([0 1]);
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0);
pbaspect([10,1,1]);
colormap(CUB); caxis([-2 2]);
colorbar('Ticks',linspace(minZx,maxZx,2),'FontSize',fSiz-5);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-0.05 1.05]);

if max(panel23_3C_Y,[],'all') > eps
    maxZy = floor(10*max(panel23_3C_Y,[],'all'))*0.1;
else
    maxZy = ceil(10*max(panel23_3C_Y,[],'all'))*0.1;
end
if min(panel23_3C_Y,[],'all') < -eps
    minZy = ceil(10*min(panel23_3C_Y,[],'all'))*0.1;
else
    minZy = floor(10*min(panel23_3C_Y,[],'all'))*0.1;
end

subplot(3,2,4)
surf(a,c23,panel23_3C_Y,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
yticks([0 1]);
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0); 
pbaspect([10,1,1]);
colormap(CUB); caxis([-2 2]);
colorbar('Ticks',linspace(minZy,maxZy,2),'FontSize',fSiz-5);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-0.05 1.05]);

if max(panel12_3C_X,[],'all') > eps
    maxZx = floor(10*max(panel12_3C_X,[],'all'))*0.1;
else
    maxZx = ceil(10*max(panel12_3C_X,[],'all'))*0.1;
end
if min(panel12_3C_X,[],'all') < -eps
    minZx = ceil(10*min(panel12_3C_X,[],'all'))*0.1;
else
    minZx = floor(10*min(panel12_3C_X,[],'all'))*0.1;
end

subplot(3,2,1)
surf(a,c12,panel12_3C_X,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
y1 = ylabel('$c$','fontsize',fSiz,...
    'HorizontalAlignment','center','Interpreter','latex');
y1.Position(1) = y1.Position(1) + val;
xlabel('\alpha','fontweight','bold','fontsize',fSiz);
yticks([-1 0]);
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0);
pbaspect([10,1,1]);
colormap(CUB); caxis([-2 2]);
colorbar('Ticks',linspace(minZx,maxZx,2),'FontSize',fSiz-5);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.05 0.05]);

if max(panel12_3C_Y,[],'all') > eps
    maxZy = floor(10*max(panel12_3C_Y,[],'all'))*0.1;
else
    maxZy = ceil(10*max(panel12_3C_Y,[],'all'))*0.1;
end
if min(panel12_3C_Y,[],'all') < -eps
    minZy = ceil(10*min(panel12_3C_Y,[],'all'))*0.1;
else
    minZy = floor(10*min(panel12_3C_Y,[],'all'))*0.1;
end

subplot(3,2,2)
surf(a,c12,panel12_3C_Y,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
yticks([-1 0]);
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0);
pbaspect([10,1,1]);
colormap(CUB); caxis([-2 2]);
colorbar('Ticks',linspace(minZy,maxZy,2),'FontSize',fSiz-5);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.05 0.05]);

if max(panel31_3C_X,[],'all') > eps
    maxZx = floor(10*max(panel31_3C_X,[],'all'))*0.1;
else
    maxZx = ceil(10*max(panel31_3C_X,[],'all'))*0.1;
end
if min(panel31_3C_X,[],'all') < -eps
    minZx = ceil(10*min(panel31_3C_X,[],'all'))*0.1;
else
    minZx = floor(10*min(panel31_3C_X,[],'all'))*0.1;
end

subplot(3,2,5)
surf(a,c13,panel31_3C_X,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
yticks([-1 1]);
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0); 
pbaspect([10,1,1]);
colormap(CUB); caxis([-2 2]);
colorbar('Ticks',linspace(minZx,maxZx,2),'FontSize',fSiz-5);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.1 1.1]);

if max(panel31_3C_Y,[],'all') > eps
    maxZy = floor(10*max(panel31_3C_Y,[],'all'))*0.1;
else
    maxZy = ceil(10*max(panel31_3C_Y,[],'all'))*0.1;
end
if min(panel31_3C_Y,[],'all') < -eps
    minZy = ceil(10*min(panel31_3C_Y,[],'all'))*0.1;
else
    minZy = floor(10*min(panel31_3C_Y,[],'all'))*0.1;
end

subplot(3,2,6)
surf(a,c13,panel31_3C_Y,'EdgeColor','none');
view(0,90);
ax = gca;
ax.FontSize = fSiz;
yticks([-1 1]);
xticks([-pi/2 0 pi/2]); xticklabels({'-\pi/2','0','\pi/2'});
xtickangle(0);
pbaspect([10,1,1]);
colormap(CUB); caxis([-2 2]);
colorbar('Ticks',linspace(minZy,maxZy,2),'FontSize',fSiz-5);
xlim([-pi/2-pi/36 pi/2+pi/36]); ylim([-1.1 1.1]);


%% GAITS FOR 3C
% Here, we just alternate between each panel to obtain a gait estimate.

clear gcirc3C_12 gcirc3C_23 gcirc3C_31;

% Load the data after switching between the respective gaits
% load('sysf_hybridC_3c1r__shchf_hybridC_Square.mat'); gcirc3C_23 = p.G_locus_full{1,1}.G(:,1:2);
% load('sysf_hybridC_3c1r__shchf_hybridC_Square.mat'); gcirc3C_12 = p.G_locus_full{1,1}.G(:,1:2);
% load('sysf_hybridC_3c1r__shchf_hybridC_Square.mat'); gcirc3C_31 = p.G_locus_full{1,1}.G(:,1:2);
% Let's save the data we have locally:
% save('contactPair122331Displacement_3hybC.mat',"gcirc3C_12","gcirc3C_23","gcirc3C_31");

% Colors:
cR = [239,59,44]/255; cB = [66,146,198]/255; cG = [65,171,93]/255;

fontS = 10;

% Load the data:
load('contactPair122331Displacement_3hybC.mat');

% Since the legs are 120degs apart in the 3 contact case, let's make the BL
% the corresponding projection into the body's x-axis:
bl = (2/3);

% Normalize the displacements:
gcirc3C_12 = gcirc3C_12; gcirc3C_23 = gcirc3C_23;
gcirc3C_31 = gcirc3C_31;

% Plot the trajectory data:
f = figure();
h1 = plot(gcirc3C_12(:,1),gcirc3C_12(:,2),'color',cG,'linewidth',2);
grid on; hold on;
yline(0,'k--','LineWidth',1.2); 
xline(0,'k--','LineWidth',1.2);
h2 = plot(gcirc3C_23(:,1),gcirc3C_23(:,2),'color',cR,'linewidth',2);
h3 = plot(gcirc3C_31(:,1),gcirc3C_31(:,2),'color',cB,'linewidth',2);
axis equal; axis([-1.15 1.7 -1.4 0.9]); 
ylabel('Displacement (BLs)','fontweight','bold','fontsize',fontS);
xlabel('Displacement (BLs)','fontweight','bold','fontsize',fontS);
title('Displacements for square gaits in each panel');
legend([h1 h2 h3], 'c12', 'c23', 'c31', 'location', 'northwest');































































































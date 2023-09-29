% This script generates the contact interpolation functions.

clear all; close all; clc;

% Colors:
% C1 = (1/255)*[239,59,44]; C2 = (1/255)*[65,171,93]; C3 = (1/255)*[66,146,198];
C1 = [215,25,28]*(1/255); C2 = [77,172,38]*(1/255); C3 = [44,123,182]*(1/255);

% Define the contact domain:
c = linspace(-1.25,1.25,200);

% Initialize the variables:
f1_2 = nan(size(c)); f2_2 = nan(size(c)); 
f1_3 = nan(size(c)); f2_3 = nan(size(c)); f3_3 = nan(size(c));

% Get the data in a loop:
for i = 1:length(c)
    % Compute the interpolation function - 2 hybC case:
    f1_2(i) = contactMap(c(i),1,2); f2_2(i) = contactMap(c(i),2,2);
    % Compute the interpolation function - 3 hybC case:
    f1_3(i) = contactMap(c(i),1,3); f2_3(i) = contactMap(c(i),2,3); 
    f3_3(i) = contactMap(c(i),3,3);
end

% Label font size:
fontS = 15;

% Make the plots and save:
f = figure(); %'WindowState','maximized'
% h1 = subplot(2,1,1);
plot(c,f1_2,'color',C1,'linewidth',3,'LineStyle','--'); % h1 = 
grid on; hold on;
ax = gca;
ax.FontSize = fontS;
plot(c,f2_2,'color',C2,'linewidth',3,'LineStyle','--'); % h2 = 
xticks([-1, 1]); yticks([0, 1]);
xline([-1,1],'k--','linewidth',1.2);
% ylabel('$f(c)$','fontsize',fontS,'Interpreter','latex','fontweight','bold'); %,'fontweight','bold'
pbaspect([3 1 1]);
axis([-1 1 0 1]);

% disp(h1.Position)

% clear f1_2 f2_2 Col2
% 
% discN = 1000;
% 
% a = linspace(-pi/2,pi/2,discN);
% c = linspace(-1,1,discN);
% for i = 1:length(c)
%     f1_2(i) = contactMap(c(i),1,2); f2_2(i) = contactMap(c(i),2,2);
% end
% 
% 
% % Make everything into a matrix:
% [A,C] = meshgrid(a,c);
% f1_2 = repmat(f1_2,discN,1); f2_2 = repmat(f2_2,discN,1);
% 
% % Get the color:
% Col2 = permute(C1,[1 3 2]).*f1_2' + permute(C2,[1 3 2]).*f2_2';
% 
% h2 = subplot(2,1,2);
% surf(C,A,f1_2 + f2_2 - 1,Col2,'EdgeColor','none');
% grid on; hold on;
% ax = gca;
% ax.FontSize = fontS;
% % % surf(a,c,f2_2,Col2,'color',C3,'linewidth',3);
% xticks([-1,1]); 
% yticks([0,1]);
% xline([-1,1],'k--','linewidth',1.2);
% pbaspect([3 0.6 1]);
% % ylabel('$\sum f(c)$','fontsize',fontS,'Interpreter','latex','fontweight','bold'); %,'fontweight','bold'
% axis([-1 1 0 1]);
% view(0,90);
% 
% % sgtitle('Contact Strength','fontweight','bold','fontsize',fontS);
% 
% shift_Left = 0.0; shift_Up = -0.08; 
% 
% disp(h2.Position)
% h2.Position = h2.Position - [shift_Left,shift_Up,0,0];

% exportgraphics(f,'ContactInterpolationFunction_sine_2hybC_Strength.png','Resolution',600);
% print('-painters','-dsvg','ContactInterpolationFunction_sine_2hybC_Strength.svg');
% 
% clear c;
% c = linspace(-1.25,1.25,200);


f = figure(); % 'WindowState','maximized'
% h3 = subplot(2,1,1);
plot(c,f1_3,'color',C1,'linewidth',3,'LineStyle','--');
grid on; hold on;
ax = gca;
ax.FontSize = fontS;
plot(c,f2_3,'color',C2,'linewidth',3,'LineStyle','--');
plot(c,f3_3,'color',C3,'linewidth',3,'LineStyle','--');
pbaspect([3 1 1]);
xticks([-1, 0, 1]); 
yticks([0, 1]);
xline([-1,0,1],'k--','linewidth',1.2);
% ylabel('$f(c)$','fontsize',fontS,'Interpreter','latex','fontweight','bold'); %,'fontweight','bold'
axis([-1 1 0 1]);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIAL PURPLE
C2 = [96,0,220]*(1/255);

% Make the plots and save:
f = figure(); %'WindowState','maximized'
% h1 = subplot(2,1,1);
plot(c,f1_2,'color',C1,'linewidth',3,'LineStyle','--'); % h1 = 
grid on; hold on;
ax = gca;
ax.FontSize = 2*fontS;
plot(c,f2_2,'color',C2,'linewidth',3,'LineStyle','--'); % h2 = 
% xticks([-1, 1]); yticks([-1, 1]);
xticks([]); yticks([]);
xline([-1,1],'k--','linewidth',1.2);
% ylabel('$f(c)$','fontsize',fontS,'Interpreter','latex','fontweight','bold'); %,'fontweight','bold'
pbaspect([1 1 1]);
axis([-1 1 0 1]);

clear c;
c = linspace(-1.25,1.25,200);


f = figure(); % 'WindowState','maximized'
% h3 = subplot(2,1,1);
plot(c,f1_3,'color',C1,'linewidth',3,'LineStyle','--');
grid on; hold on;
ax = gca;
ax.FontSize = fontS;
plot(c,f2_3,'color',C2,'linewidth',3,'LineStyle','--');
plot(c,f3_3,'color',C3,'linewidth',3,'LineStyle','--');
pbaspect([3 1 1]);
xticks([-1, 0, 1]); 
yticks([0, 1]);
xline([-1,0,1],'k--','linewidth',1.2);
% ylabel('$f(c)$','fontsize',fontS,'Interpreter','latex','fontweight','bold'); %,'fontweight','bold'
axis([-1 1 0 1]);

% disp(h3.Position)

% clear f1_3 f2_3 f3_3 Col3
% 
% C1 = [215,25,28]*(1/255); C2 = [77,172,38]*(1/255); C3 = [44,123,182]*(1/255);
% 
% discN = 1000;
% 
% a = linspace(-pi/2,pi/2,discN);
% c = linspace(-1,1,discN);
% for i = 1:length(c)
%     f1_3(i) = contactMap(c(i),1,3); f2_3(i) = contactMap(c(i),2,3); 
%     f3_3(i) = contactMap(c(i),3,3);
% end
% 
% 
% % Make everything into a matrix:
% [A,C] = meshgrid(a,c);
% f1_3 = repmat(f1_3,discN,1); f2_3 = repmat(f2_3,discN,1); f3_3 = repmat(f3_3,discN,1);
% 
% % Get the color:
% Col3 = permute(C1,[1 3 2]).*f1_3' + permute(C2,[1 3 2]).*f2_3' +...
%     permute(C3,[1 3 2]).*f3_3';
% 
% h4 = subplot(2,1,2);
% surf(C,A,f1_3 + f2_3 + f3_3 - 1,Col3,'EdgeColor','none');
% grid on; hold on;
% ax = gca;
% ax.FontSize = fontS;
% pbaspect([3 0.6 1]);
% xticks([]); 
% yticks([0, 1]);
% xline([-1,0,1],'k--','linewidth',1.2);
% axis([-1 1 0 1]);
% % ylabel('\bf{$\sum f(c)$}','fontsize',fontS,'Interpreter','latex'); %,'fontweight','bold'
% view(0,90);
% 
% disp(h4.Position)
% h4.Position = h4.Position - [shift_Left,shift_Up,0,0];

% exportgraphics(f,'ContactInterpolationFunction_sine_3hybC_Strength.png','Resolution',600);
% print('-painters','-dsvg','ContactInterpolationFunction_sine_3hybC_Strength.svg');

% sgtitle('Contact Strength','fontweight','bold','fontsize',fontS);

% %%
% 
% clear f1_2 f2_2 f1_3 f2_3 f3_3 Col2 Col3
% 
% C1 = [215,25,28]*(1/255); C2 = [77,172,38]*(1/255); C3 = [44,123,182]*(1/255);
% 
% discN = 1000;
% 
% a = linspace(-pi/2,pi/2,discN);
% c = linspace(-1,1,discN);
% 
% % Get the data in a loop:
% for i = 1:length(c)
%     % Compute the interpolation function - 2 hybC case:
%     f1_2(i) = contactMap(c(i),1,2); f2_2(i) = contactMap(c(i),2,2);
%     % Compute the interpolation function - 3 hybC case:
%     f1_3(i) = contactMap(c(i),1,3); f2_3(i) = contactMap(c(i),2,3); 
%     f3_3(i) = contactMap(c(i),3,3);
% end
% 
% 
% % Make everything into a matrix:
% [A,C] = meshgrid(a,c);
% f1_2 = repmat(f1_2,discN,1); f2_2 = repmat(f2_2,discN,1);
% f1_3 = repmat(f1_3,discN,1); f2_3 = repmat(f2_3,discN,1); f3_3 = repmat(f3_3,discN,1);
% 
% % Get the color:
% Col2 = permute(C1,[1 3 2]).*f1_2' + permute(C2,[1 3 2]).*f2_2';
% Col3 = permute(C1,[1 3 2]).*f1_3' + permute(C2,[1 3 2]).*f2_3' +...
%     permute(C3,[1 3 2]).*f3_3';
% 
% %%
% fontS = 15;
% 
% f = figure();
% surf(A,C,f1_2 + f2_2 - 1,Col2,'EdgeColor','none');
% grid on; hold on;
% ax = gca;
% ax.FontSize = fontS;
% % % surf(a,c,f2_2,Col2,'color',C3,'linewidth',3);
% yticks([-1, 1]); xticks([0, 1]);
% yline([-1,1],'k--','linewidth',1.2);
% pbaspect([1 5 1]);
% % title('Contact Strength','fontweight','bold','fontsize',fontS-5);
% axis([0 1 -1 1]);
% view(0,90);
% % exportgraphics(f,'ContactInterpolationFunction_sine_2hybC_Strength.png','Resolution',500);
% 
% f = figure();
% surf(A,C,f1_3 + f2_3 + f3_3 - 1,Col3,'EdgeColor','none');
% grid on; hold on;
% ax = gca;
% ax.FontSize = fontS;
% pbaspect([1 5 1]);
% yticks([-1, 0, 1]); xticks([0, 1]);
% yline([-1,0,1],'k--','linewidth',1.2);
% % title('Contact Strenght','fontweight','bold','fontsize',fontS);
% axis([0 1 -1 1]);
% view(0,90);
% % exportgraphics(f,'ContactInterpolationFunction_sine_3hybC_Strength.png','Resolution',500);

%% 2C BASELINE GAIT FOR OPTIMIZATION:

% Define your gait phase:
phi = linspace(0,1,1001);

% Initialize your variables:
C = nan(size(phi));
alpha = nan(size(phi));

% Define the contact gait:
C(phi <= 0.25) = -1;
C(phi > 0.25 & phi <= 0.5) = 8*(phi(phi > 0.25 & phi <= 0.5) - 0.25) - 1;
C(phi > 0.5 & phi <= 0.75) = 1;
C(phi > 0.75 & phi <= 1) = -8*(phi(phi > 0.75 & phi <= 1) - 0.75) + 1;

% Define the swing gait:
alpha(phi >= 0.75) = C(phi <= 0.25);
alpha(phi < 0.75) = C(phi > 0.25);

% Define your plotting parameters here:
fSiz = 90; linW = 10;

% Plot the data:
f = figure('WindowState','maximized'); % 'WindowState','maximized'
% yyaxis left
plot(phi(phi <= 0.25),alpha(phi <= 0.25),...
    'Color',C1,'LineWidth',linW,'LineStyle','-');
hold on;
plot(phi(phi > 0.25 & phi <= 0.5),alpha(phi > 0.25 & phi <= 0.5),...
    'Color','k','LineWidth',linW,'LineStyle','--');
plot(phi(phi > 0.5 & phi <= 0.75),alpha(phi > 0.5 & phi <= 0.75),...
    'Color',C2,'LineWidth',linW,'LineStyle','-');
plot(phi(phi > 0.75 & phi <= 1),alpha(phi > 0.75 & phi <= 1),...
    'Color','k','LineWidth',linW,'LineStyle','--');
    xlabel('Time Period \tau','FontSize',fSiz);
ylabel('\alpha','FontSize',fSiz);
yticks([-1,1]); yticklabels({'\alpha^-','\alpha^+'});
xticks([0,0.25,0.5,0.75,1]);
% xline([0,0.25,0.5,0.75,1],'k--','LineWidth',1.2);
ylim([-1.1 1.1]);

pbaspect([1 0.2 1]);

text(0.135,0,'c = -1','Color',C1,'FontSize',fSiz); 
text(0.640,0,'c = 1','Color',C2,'FontSize',fSiz);

% yyaxis right
% plot(phi,C,'k','LineWidth',linW);
% ylabel('Contact State','FontSize',fSiz,'Color','k');
% yticks([-1,1]); yticklabels({'\{1,0\}','\{0,1\}'});
% ylim([-1.1 1.1]);

ax = gca;
ax.FontSize = fSiz;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';

% exportgraphics(f,'Gait_2C_Opt_newC_V2_Overleaf.png','Resolution',500);

%%

% V1 function
function fi = contactMap(c,i,n)



    % In V1, first and last contact states are formed by 1 and n
    % respectively.

    % Get the spacing (V1):
    s = 2/(n-1);

    % Get the sinusoidal term coefficient:
    coeff = pi/s;
    
    % Switch the case based on which interpolation function is needed.
    switch 1


        case i == 1 % first contact state
            
            if c <= -1
                fi = 1;
            elseif c <= -1+s && c > -1
                fi = ( 1 + cos(coeff*(c+1)) )*0.5;
            else
                fi = 0;
            end

        case i == n % last contact state

            if c <= -1 + (n-2)*s
                fi = 0;
            elseif c <= 1 && c > -1 + (n-2)*s
                fi = ( 1 - cos(coeff*(c + 1 - (n-2)*s)) )*0.5;
            else
                fi = 1;
            end

        case i > 1 && i < n % intermediate contact states

            if c <= -1 + (i-2)*s
                fi = 0;
            elseif c <= -1 + i*s && c > -1 + (i-2)*s
                fi = ( 1 - cos(coeff*(c + 1 - (i-2)*s)) )*0.5;
            else
                fi = 0;
            end


    end



end

% %%
% 
% %%% interpolation function - version 2 (contact connects on itself)
% 
% function fi = contactMap(c,i,n)
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % In V2, first and last contact states are formed by 1. All other
%     % contact states are intermediate.
% 
%     % Get the spacing (V2):
%     s = 2/n;
% 
%     % Get the sinusoidal term coefficient:
%     coeff = pi/s;
%     
%     % Switch the case based on which interpolation function is needed.
%     switch 1
% 
% 
%         case i == 1 % first contact state
%             
%             if c <= -1
%                 fi = 1;
%             elseif c <= -1+s && c > -1
%                 fi = ( 1 + cos(coeff*(c+1)) )*0.5;
%             elseif c > -1+s && c <= -1 + (n-1)*s
%                 fi = 0;
%             elseif c > -1 + (n-1)*s && c <= 1
%                 fi = ( 1 - cos(coeff*(c + 1 -(n-1)*s)) )*0.5;
%             else % greater than or equal to 1 case
%                 fi = 1;
%             end
% 
%         case i ~= 1 % all other contact states
% 
%             if c <= -1 + (i-2)*s
%                 fi = 0;
%             elseif c <= -1 + i*s && c > -1 + (i-2)*s
%                 fi = ( 1 - cos(coeff*(c + 1 - (i-2)*s)) )*0.5;
%             elseif c > -1 + i*s
%                 fi = 0;
%             end
% 
% 
%     end
% 
% 
% 
% end


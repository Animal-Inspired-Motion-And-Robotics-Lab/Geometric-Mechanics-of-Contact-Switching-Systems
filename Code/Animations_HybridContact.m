%% GENERAL PARAMETERS

clear all; close all; clc;

% Color map:
col = (1/255)*[215,25,28;
            96,0,220;
            44,123,182;
            77,175,74;
            255,127,0];

%% COMPLEX CASE: 4-leg system with mid-start gait image

% For now, we can use standard ankle limits of 90 degs on each side -- fix
% this later.

% Initialize the simulation structure and parameters:
sim = [];
numL = 4; % number of legs
sim.phi_pin_ord = [1 4 2 3 1]; % each submanifold visited in chronological order from left to right
sim.phi_range = [-45 0 90 0 -90; % the rotor ranges for each submanifold
            0 90 0 -90 -45]; sim.phi_range = deg2rad(sim.phi_range);
sim.T = [1 0.5]; % all legs have the same swinging time (left) and switching time (right)
% The first encodes the time needed to swing between the rotor limits.
sim.off = (1/numL)*2*pi*(linspace(1,numL,numL) - 1); % symmetrically distributed legs

% The final time should be computed by figuring out the time taken for
% swinging and switching
alpha_dot = pi/sim.T(1); % the provided time to swing is for the full rotor range
                     % convert this into the constant swing rate
T_alpha = abs(diff(sim.phi_range,1,1))/(alpha_dot);
T_beta = sim.T(2)*ones(size(T_alpha)); % the contact time should be the same
if sim.phi_pin_ord(end) == sim.phi_pin_ord(1)
    T_leg = T_alpha + [T_beta(1:end-1), 0];
else
    T_leg = T_alpha + T_beta;
end
% sim.phi_tf = 3*sum(T_leg); % Let's simulate it for 3 cycles
sim.phi_tf = 1*sum(T_leg); % Let's simulate it for 1 cycle

% Set a boundary dimension (remember the leg is 1 unit long)
legL = 0.01; % 0.25

% Initialize the video structure and parameters:
vid = [];
vid.frame_r = 100;
vid.name = "four_leg_split_1423";
vid.shch_name = "four_leg_split_1423_shape";

% Get the computations
traj = ComputeSingleDisk_v2(sim,vid);

% Add a leg circle scattering size:
circS = 250;

% Add this to a separate time structure for processing in animate function
tS = [];
tS.T_leg = T_leg;
tS.phi_pin_ord = sim.phi_pin_ord;

Xpos = traj(5:4+numL,1:end-1); Ypos = traj(4+numL+1:4+2*numL,1:end-1);

limX = [min(Xpos,[],'all') - legL, max(Xpos,[],'all') + legL];
limY = [min(Ypos,[],'all') - legL, max(Ypos,[],'all') + legL];

% Trajectory color:
traj_c = zeros(1,3); % black

% Provide the leg color:
if size(col,1) < numL
    leg_c = interp1([0,1],[col(1,:); col(end,:)],linspace(0,1,numL));
else
    leg_c = col(1:numL,:);
end

% Animate the system -- also animate the shape space
AnimateSingleDisk(traj,traj_c,leg_c,vid,limX,limY,circS,tS,sim.off); % ,tS

%% SIMPLE CASE: 2-leg system with different swing and switch time-periods

%%%%% 1
sim = [];
numL = 2;
sim.phi_pin_ord = [1 2];

amp = 0.72;
sim.phi_range = [-amp, amp;
                  amp, -amp];

sim.T = [0.45838 0.1];
sim.off = (1/numL)*2*pi*(linspace(1,numL,numL) - 1);
alpha_dot = pi/sim.T(1);
T_alpha = abs(diff(sim.phi_range,1,1))/(alpha_dot);
T_beta = sim.T(2)*ones(size(T_alpha));
if sim.phi_pin_ord(end) == sim.phi_pin_ord(1)
    T_leg = T_alpha + [T_beta(1:end-1), 0];
else
    T_leg = T_alpha + T_beta;
end
sim.phi_tf = 4; % manually chosen

legL = 0.01;

vid = [];
vid.frame_r = 100;
vid.name = ['two_leg_normal_12_optamp', num2str(amp), '_Fig4'];
vid.shch_name = ['two_leg_normal_12_optamp', num2str(amp), '_Fig4_shape'];

traj = ComputeSingleDisk_v2(sim,vid);

circS = 250;

tS = [];
tS.T_leg = T_leg;
tS.phi_pin_ord = sim.phi_pin_ord;

Xpos = traj(5:4+numL,1:end-1); Ypos = traj(4+numL+1:4+2*numL,1:end-1);
limX = [min(Xpos,[],'all') - legL, max(Xpos,[],'all') + legL];
limY = [min(Ypos,[],'all') - legL, max(Ypos,[],'all') + legL];

traj_c = zeros(1,3);

leg_c = col(1:2,:);

AnimateSingleDisk(traj,traj_c,leg_c,vid,limX,limY,circS,tS,sim.off); % ,tS

%%%%% 2

amp = 1.23;
sim.phi_range = [-amp, amp;
                  amp, -amp];

sim.T = [0.7815 1];
alpha_dot = pi/sim.T(1);
T_alpha = abs(diff(sim.phi_range,1,1))/(alpha_dot);
T_beta = sim.T(2)*ones(size(T_alpha));
if sim.phi_pin_ord(end) == sim.phi_pin_ord(1)
    T_leg = T_alpha + [T_beta(1:end-1), 0];
else
    T_leg = T_alpha + T_beta;
end

vid.name = ['two_leg_normal_12_optamp', num2str(amp), '_Fig4'];
vid.shch_name = ['two_leg_normal_12_optamp', num2str(amp), '_Fig4_shape'];

traj = ComputeSingleDisk_v2(sim,vid);

circS = 250;

tS = [];
tS.T_leg = T_leg;
tS.phi_pin_ord = sim.phi_pin_ord;

Xpos = traj(5:4+numL,1:end-1); Ypos = traj(4+numL+1:4+2*numL,1:end-1);
limX = [min(Xpos,[],'all') - legL, max(Xpos,[],'all') + legL];
% limY = [min(Ypos,[],'all') - legL, max(Ypos,[],'all') + legL];
% The y-limit will be used without change as before to make the comparison
% easier.

AnimateSingleDisk(traj,traj_c,leg_c,vid,limX,limY,circS,tS,sim.off); % ,tS

%%

%%%%% EXAMPLE CASE

% number of gait cycles
numC = 2;

numL = 2;
amp = pi/4;
sim.phi_pin_ord = [1 2];
sim.phi_range = [-amp, amp;
                  amp, -amp];
sim.T = [1 1]; % computed these switching times from the paper optimization section
sim.off = (1/numL)*2*pi*(linspace(1,numL,numL) - 1);

alpha_dot = pi/sim.T(1);
T_alpha = abs(diff(sim.phi_range,1,1))/(alpha_dot);
T_beta = sim.T(2)*ones(size(T_alpha));
if sim.phi_pin_ord(end) == sim.phi_pin_ord(1)
    T_leg = T_alpha + [T_beta(1:end-1), 0];
else
    T_leg = T_alpha + T_beta;
end
sim.phi_tf = numC*sum(T_leg);

legL = 0.01;

vid.frame_r = 100;
vid.name = "two_leg_normal_12_Fig3_example";
vid.shch_name = "two_leg_normal_12_Fig3_example_shape";

traj = ComputeSingleDisk_v2(sim,vid);

traj_c = zeros(1,3);
leg_c = col(1:2,:);

circS = 250;

tS = [];
tS.T_leg = T_leg;
tS.phi_pin_ord = sim.phi_pin_ord;

Xpos = traj(5:4+numL,1:end-1); Ypos = traj(4+numL+1:4+2*numL,1:end-1);

limX = [min(Xpos,[],'all') - legL, max(Xpos,[],'all') + legL];
limY = [min(Ypos,[],'all') - legL, max(Ypos,[],'all') + legL];

AnimateSingleDisk(traj,traj_c,leg_c,vid,limX,limY,circS,tS,sim.off);

%%

%%%%% INFINITESIMAL STRENGTH CASE

% number of gait cycles
numC = 1;

numL = 2;

%%%%%% point 1 at alpha = -45degs
pt = -pi/2;
amp = pi/24; % six time smaller amplitude than the example
sim.phi_pin_ord = [1 2];
sim.phi_range = [pt-amp, pt+amp;
                 pt+amp, pt-amp];
sim.T = [1 1]; % computed these switching times from the paper optimization section
sim.off = (1/numL)*2*pi*(linspace(1,numL,numL) - 1);

alpha_dot = pi/sim.T(1);
T_alpha = abs(diff(sim.phi_range,1,1))/(alpha_dot);
T_beta = sim.T(2)*ones(size(T_alpha));
if sim.phi_pin_ord(end) == sim.phi_pin_ord(1)
    T_leg = T_alpha + [T_beta(1:end-1), 0];
else
    T_leg = T_alpha + T_beta;
end
sim.phi_tf = numC*sum(T_leg);

legL = 0.01;

vid.frame_r = 100;
vid.name = "two_leg_inf1_12";
vid.shch_name = "two_leg_inf1_12_shape";

traj = ComputeSingleDisk_v2(sim,vid);

traj_c = zeros(1,3);
leg_c = col(1:2,:);

circS = 250;

tS = [];
tS.T_leg = T_leg;
tS.phi_pin_ord = sim.phi_pin_ord;

limOff = 1;
Xpos = traj(5:4+numL,1:end-1); Ypos = traj(4+numL+1:4+2*numL,1:end-1);
limX = [min(Xpos,[],'all') - legL - limOff, max(Xpos,[],'all') + legL + limOff];
limY = [min(Ypos,[],'all') - legL - limOff, max(Ypos,[],'all') + legL + limOff];

AnimateSingleDisk(traj,traj_c,leg_c,vid,limX,limY,circS,tS,sim.off);

%%%%%%%%%%%%% point 2 at alpha = 22.5degs
pt = -pi/2*(0.5);
amp = pi/24; % six time smaller amplitude than the example
sim.phi_pin_ord = [1 2];
sim.phi_range = [pt-amp, pt+amp;
                 pt+amp, pt-amp];
sim.T = [1 1]; % computed these switching times from the paper optimization section
sim.off = (1/numL)*2*pi*(linspace(1,numL,numL) - 1);

alpha_dot = pi/sim.T(1);
T_alpha = abs(diff(sim.phi_range,1,1))/(alpha_dot);
T_beta = sim.T(2)*ones(size(T_alpha));
if sim.phi_pin_ord(end) == sim.phi_pin_ord(1)
    T_leg = T_alpha + [T_beta(1:end-1), 0];
else
    T_leg = T_alpha + T_beta;
end
sim.phi_tf = numC*sum(T_leg);

legL = 0.01;

vid.frame_r = 100;
vid.name = "two_leg_inf2_12";
vid.shch_name = "two_leg_inf2_12_shape";

traj = ComputeSingleDisk_v2(sim,vid);

traj_c = zeros(1,3);
leg_c = col(1:2,:);

circS = 250;

tS = [];
tS.T_leg = T_leg;
tS.phi_pin_ord = sim.phi_pin_ord;

limOff = 1;
Xpos = traj(5:4+numL,1:end-1); Ypos = traj(4+numL+1:4+2*numL,1:end-1);
limX = [min(Xpos,[],'all') - legL - limOff, max(Xpos,[],'all') + legL + limOff];
limY = [min(Ypos,[],'all') - legL - limOff, max(Ypos,[],'all') + legL + limOff];

AnimateSingleDisk(traj,traj_c,leg_c,vid,limX,limY,circS,tS,sim.off);

%%%%%%%%%%%%% point 3 at alpha = 0degs
pt = -pi/2*(0.0);
amp = pi/24; % six time smaller amplitude than the example
sim.phi_pin_ord = [1 2];
sim.phi_range = [pt-amp, pt+amp;
                 pt+amp, pt-amp];
sim.T = [1 1]; % computed these switching times from the paper optimization section
sim.off = (1/numL)*2*pi*(linspace(1,numL,numL) - 1);

alpha_dot = pi/sim.T(1);
T_alpha = abs(diff(sim.phi_range,1,1))/(alpha_dot);
T_beta = sim.T(2)*ones(size(T_alpha));
if sim.phi_pin_ord(end) == sim.phi_pin_ord(1)
    T_leg = T_alpha + [T_beta(1:end-1), 0];
else
    T_leg = T_alpha + T_beta;
end
sim.phi_tf = numC*sum(T_leg);

legL = 0.01;

vid.frame_r = 100;
vid.name = "two_leg_inf3_12";
vid.shch_name = "two_leg_inf3_12_shape";

traj = ComputeSingleDisk_v2(sim,vid);

traj_c = zeros(1,3);
leg_c = col(1:2,:);

circS = 250;

tS = [];
tS.T_leg = T_leg;
tS.phi_pin_ord = sim.phi_pin_ord;

limOff = 1;
Xpos = traj(5:4+numL,1:end-1); Ypos = traj(4+numL+1:4+2*numL,1:end-1);
limX = [min(Xpos,[],'all') - legL - limOff, max(Xpos,[],'all') + legL + limOff];
limY = [min(Ypos,[],'all') - legL - limOff, max(Ypos,[],'all') + legL + limOff];

AnimateSingleDisk(traj,traj_c,leg_c,vid,limX,limY,circS,tS,sim.off);

%% COMPLEX CASE: 3-leg system (Paper Fig 7)
sim = [];
numL = 3;
sim.phi_pin_ord = [1 3 2 3 1];
sim.phi_range = [0 45 90 -45 -90;
                45 90 -45 -90 0]; 
sim.phi_range = deg2rad(sim.phi_range);
sim.T = [1 0.5];
sim.off = (1/numL)*2*pi*(linspace(1,numL,numL) - 1);

alpha_dot = pi/sim.T(1);
T_alpha = abs(diff(sim.phi_range,1,1))/(alpha_dot);
T_beta = sim.T(2)*ones(size(T_alpha));
if sim.phi_pin_ord(end) == sim.phi_pin_ord(1)
    T_leg = T_alpha + [T_beta(1:end-1), 0];
else
    T_leg = T_alpha + T_beta;
end
sim.phi_tf = 1*sum(T_leg);

legL = 0.1; % 0.25

vid = [];
vid.frame_r = 100;
vid.name = "three_leg_split_13231_Fig7";
vid.shch_name = "three_leg_split_13231_Fig7_shape";

circS = 250;

traj = ComputeSingleDisk_v2(sim,vid);

tS = [];
tS.T_leg = T_leg;
tS.phi_pin_ord = sim.phi_pin_ord;

Xpos = traj(5:4+numL,1:end-1); Ypos = traj(4+numL+1:4+2*numL,1:end-1);
limX = [min(Xpos,[],'all') - legL, max(Xpos,[],'all') + legL];
limY = [min(Ypos,[],'all') - legL, max(Ypos,[],'all') + legL];

traj_c = zeros(1,3);

% leg_c = (1/255)*[215, 25, 28;
%                 96, 0 ,220;
%                 44, 123, 182];
leg_c = col(1:3,:);

AnimateSingleDisk(traj,traj_c,leg_c,vid,limX,limY,circS,tS,sim.off); %,circS,tS

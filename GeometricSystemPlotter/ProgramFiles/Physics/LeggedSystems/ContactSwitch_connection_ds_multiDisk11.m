% Calculate the local connection for a 2-legged multi-disk system
function [A,h,J,J_full,omega] = ContactSwitch_connection_ds_multiDisk11(geometry,~,jointangles)

    % let's set the link length as 1
    l = 1;

    % Get Pfaffian constraint relating limb velocities, sprawl-velocity, 
    % and left leg's rotational velocity
    Omega = @(a1,a2,g,l) l*[-(sin(a1-a2+2*g)+cos(a1)-cos(a1+2*g)), 2*cos(a1), 0, -(cos(a1)-cos(a1+2*g));
                            -(cos(a1-a2+2*g)-sin(a1)+sin(a1+2*g)+1), -2*(sin(a1)-1), 1, -(sin(a1+2*g)-sin(a1)+1)];

    % Now, get the connection for the fully controlled case (gamma is also
    % controlled)
    A_full = @(g) [0, 0, 2*l*cos(g);
                    -l*sin(g), -l*sin(g), 0;
                    1, 1, 0];

    % Now, we compute the pfaffian constraint matrix and A matrix for the
    % given joint angles:
    omega = Omega(jointangles(1),jointangles(2),jointangles(3),l);
    A = A_full(jointangles(3));

    % Return empty vectors/matrices for h, J and J full
    h = zeros(3,3);
    J = zeros(3,3); J_full = zeros(3,6);

    %     % Get the map from limb rotations to sprawl and foot-rotational
    %     % velocities -- NEED THIS IN THE GAIT GENERATOR
    %     sigma = @(a1,a2,g) sin(2*g) - sin(a1-a2+2*g) + cos(a1+2*g) + cos(a2-2*g);
    %     Omega_1p = @(a1,a2,g) [cos(a1), -(sin(2*g)+cos(a1+2*g));
    %                            0.5*(sin(a1-a2+2*g)+cos(a1)-cos(a1+2*g)), -0.5*(sin(a1-a2+2*g)+cos(a2)-cos(a2-2*g))]*(1/sigma(a1,a2,g));
end

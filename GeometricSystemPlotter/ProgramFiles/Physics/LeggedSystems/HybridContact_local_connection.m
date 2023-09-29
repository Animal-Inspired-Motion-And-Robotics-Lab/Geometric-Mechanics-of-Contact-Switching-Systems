function [A, h, J,J_full,Omega] = HybridContact_local_connection(geometry,physics,shapeparams)
% Calculate the local connection for a set of curvature bases
%
% Inputs:
% geometry: structure defining system geometry
%      geometry.type: how the system geometry is defined 
%         (e.g., links, curvature basis, general curvature)
%      geometry.function: map from shape variables to local backbone
%         deformation (e.g., curvature or joint angles)
%      geometry.length: total length of swimmer
% physics: structure defining system physics
%      drag_ratio: ratio of lateral to longitudinal drag
%      drag_coefficient: drag per unit length
% shapeparams: value of shape variables




% Identify what kind of system is being calculated, and use this to specify 
% how the local connection should be generated
switch geometry.type
    
    case {'n-link chain','branched chain'}
        physics_function = @HybridContact_connection_discrete;
    
    case 'multi_disk_11'
        physics_function = @ContactSwitch_connection_ds_multiDisk11;
    
end

% Call the physics function identified for the system
[A, h, J,J_full,Omega] = physics_function(geometry,physics,shapeparams);

end
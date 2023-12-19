function w_deg = convert3D_to_angV(pos0_x, pos0_z, eyePos_x, eyePos_z, vel_x, vel_z)
%%
% Based on Rokers et al., 2018 derivations.

%% example inputs. Unit is arbitrary but must be consistent.
% pos0_x = 100;
% pos0_z = 100;
% 
% vel_x = 0;
% vel_z = 135;
% 
% eyePos_x = 0;
% eyePos_z = 0;


% get intial angle of object to eye
%B = atan((pos0_z/(pos0_x-eyePos_z)));

%% main part of function

% get initial distance of object to eye
h0 = sqrt((pos0_x-eyePos_x)^2 + (pos0_z-eyePos_z)^2);

% derivative of angle change by time
dB_dt = (1/h0^2) * (vel_z*(pos0_x-eyePos_x) - vel_x*(pos0_z-eyePos_z));

w_deg = rad2deg(dB_dt); % convert to degrees
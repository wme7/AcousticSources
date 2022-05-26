function [x_i,y_i,z_i,theta_i] = build_circularArrayOfSensors(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Build cicular array of sensors 
%   Coded by Manuel A. Diaz @ Univ-Poitiers | Pprime 24.05.2022
%
%   Usage: 
%   [xS,yS,zS] = build_circularArrayOfSensors(x0,y0,z0,Radius,nSensors)
%   [xS,yS,zS] = build_circularArrayOfSensors(x0,y0,z0,Radius,nSensors,yaw,pitch,roll)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch nargin
        case 5
            xo = varargin{1};
            yo = varargin{2};
            zo = varargin{3};
            radius = varargin{4};
            nSensors =  varargin{5};
            yaw  = 0;
            pitch= 0;
            roll = 0;
        case 8
        otherwise
            disp('Usage: ')
            disp('[xS,yS,zS] = build_circularArrayOfSensors(x0,y0,z0,Radius,nSensors)')
            disp('[xS,yS,zS] = build_circularArrayOfSensors(x0,y0,z0,Radius,nSensors,yaw,pitch,roll)')
    end
    theta_i = 0:2*pi/nSensors:2*pi;
    % Sensors in a reference plane
    x_i = radius * cos(theta_i);
    y_i = radius * sin(theta_i);
    z_i = zeros( size(theta_i));
    % Rotation matrices
    T_yaw  = [1 0 0; 0 cos( yaw ) -sin( yaw ); 0 sin( yaw ) cos( yaw )];
    T_pitch= [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)];
    T_roll = [cos(roll ) -sin(roll ) 0; sin(roll ) cos(roll ) 0; 0 0 1];
    % Output final coordinates of sensors
    X = T_yaw * T_pitch * T_roll * [x_i;y_i;z_i];
    x_i = xo + X(1,:);
    y_i = yo + X(2,:);
    z_i = zo + X(3,:);
end
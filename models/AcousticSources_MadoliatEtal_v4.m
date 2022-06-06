%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Monopole, Dipole and Quadrupoles Generator v4.0
%
%          by Manuel A. Diaz @ Univ-Poitiers | Pprime - 24.05.2022
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on:
% [1] Howe, Michael S. Theory of vortex sound. No.33. 
%     Cambridge university press, 2003. 
% [2] Madoliat, Reza, Nowrouz Mohammad Nouri, and Ali Rahrovi.
%     "Equalization of acoustic source using multi-pole sources and source
%     strength estimation using inverse method." Applied Acoustics 113
%     (2016): 210-220.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

% Controling Parameters
source_type = 4;
distance = 5;

density = 1.0;
sound_speed = 1.0;
frequency = 5;
strength = 1.0;

hdistance = 5;
vdistance = 5;

rho = density;
c = sound_speed;
f = frequency;
Q = strength;
w = 2*pi*f;
k = w/c;
theta= 0:0.01:2*pi; 
%========================================================================== 
% Source Type : Monopole 
%========================================================================== 
if source_type == 1
    figure()
    for l= 1:1: length(distance)
        r= distance(l);
        for m= 1:1: length(theta)
            pressure(l,m)= abs(Q*((1i*k*rho*c)/(4*pi*r))); end
        polarplot(theta, pressure(l,:))
        hold on
    end
    hold off
end
%========================================================================== 
% Source Type : Dipole 
%========================================================================== 
if source_type == 2
    d   =  hdistance;
    figure()
    for l= 1:1: length(distance)
        r= distance(l);
        for m= 1:1: length(theta)
            f_theta= theta(m);
            pressure(l,m)= abs(((-1i*Q*rho*c*(k^2)*d)/(4*pi*r))*...
                cos(f_theta));
        end
        polarplot(theta, pressure(l,:))
        hold on
    end
    hold off
end
%==========================================================================
% Source Type : Quadrupole 
%========================================================================== 
if source_type == 3
    d= hdistance;
    D= vdistance;
    figure()
    for l= 1:1: length(distance)
        r= distance(l);
        for m= 1:1:length(theta)
            f_theta= theta(m);
            pressure(l,m)= abs((Q*rho*c*k)*(pi*r)*(k^2*d)*D*...
                cos(f_theta)*sin(f_theta));
        end
        polarplot(theta, pressure(l,:))
        hold on
    end
    hold off
end
%========================================================================== 
% Source Type : Longitudial Quadrupole 
%========================================================================== 
if source_type == 4
    d= hdistance;
    D= vdistance;
    figure()
    for l= 1:1: length(distance)
        r= distance(l);
        for m= 1:1: length(theta)
            f_theta= theta(m);
            pressure(l,m)= abs(((Q*rho*c*k)/(pi*r))*(k^2)*d*D*((cos(f_theta)^2)));
        end
        polarplot(theta, pressure(l,:))
        hold on
    end
    hold off
end
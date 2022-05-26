%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Monopole, Dipole and Quadrupoles Generator v1.0
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

% Load 
addpath ../third_party/ % third-party tools
addpath ../functions/ % sensors and ploting functions

% Media Properties
c0 = 1.0;

% Source location
xo = 0.0;
yo = 0.0;
zo = 0.0;

% Time range
t0 = 0.0;
tEnd = 1.0;

% Monopolar Source 
% Ref [1], Eq. 1.7.2
Q = @(t) 1.0*sin(2*pi*t);
psi_monopole = @(r,t) -Q(t-r/c0)./(4*pi*r);

% Dipolar Source (the x_1-directional dipole)
% Ref [1], Eq. 1.7.4 and Example 1.9.3
% NOTE: I have changed every d/dx_j by 1/c0^2*d/dt and used x_1 = r*cos(theta)
dFdt=@(t) 1.0*cos(2*pi*t); % becuase F := 1.0*sin(2*pi*t);
psi_dipole = @(r,theta,t) dFdt(t-r/c0) .* (cos(theta))./(4*pi*r*c0);

% Lateral Quadrupole Source (the (x_1,x_2)-directional quadrupole)
% Ref [1], Example 1.9.3
dT2dt2=@(t) -1.0*sin(2*pi*t); % because T := 1.0*sin(2*pi*t);
psi_quadrupole = @(r,theta,phi,t) dT2dt2(t-r/c0) .* (sin(2*theta).*cos(phi))./(8*pi*r*c0^2); 

% Numerical mesh
xa =-5; xb = 5;
ya =-4; yb = 4;
za =-3; zb = 3;
[x,y,z]=meshgrid(xa:0.1:xb,ya:0.1:yb,za:0.1:zb);
[ny,nx,nz] = size(x); 
zc = round(nz/2);

% Auxiliary variables
r=sqrt((x-xo).^2+(y-yo).^2+(z-zo).^2);
theta = atan2( (y-yo),(x-xo) ); % azimutal;
phi = atan2( (z-zo),sqrt((x-xo).^2+(y-yo).^2) ); % elevation;

% Compute acoustic pressure fields
p_monopole = psi_monopole(r,tEnd);
p_dipole = psi_dipole(r,theta,tEnd);
p_quadrupole = psi_quadrupole(r,theta,phi,tEnd);

%% Build sensors around source
[xs,ys,zs,Ts] = build_circularArrayOfSensors(xo,yo,zo,3.5,100);

% Auxiliary variables for sensors
r_s=sqrt((xs-xo).^2+(ys-yo).^2+(zs-zo).^2);
theta_s = atan2( (ys-yo),(xs-xo) ); % azimulat;
phi_s = atan2( (zs-zo),sqrt((xs-xo).^2+(ys-yo).^2) ); % elevation;

% Sample the pressure field from time t = [0,5]
t = t0:0.1:tEnd; Nt=numel(t); Ns=numel(xs);
ps_monopole   =zeros(Nt,Ns);
ps_dipole     =zeros(Nt,Ns);
ps_quadrupole =zeros(Nt,Ns);
for i = 1:Nt
    ps_monopole(i,:)   = psi_monopole(r_s,t(i));
    ps_dipole(i,:)     = psi_dipole(r_s,theta_s,t(i));
    ps_quadrupole(i,:) = psi_quadrupole(r_s,theta_s,phi_s,t(i));
end

%% Visualization 2d
figure(1); colormap whitejet; plot_source2d(x,y,z, p_monopole ,xs,ys,zs,Ts, ps_monopole , 'Monopole' );
figure(2); colormap whitejet; plot_source2d(x,y,z,  p_dipole  ,xs,ys,zs,Ts,  ps_dipole  ,  'Dipole'  );
figure(3); colormap whitejet; plot_source2d(x,y,z,p_quadrupole,xs,ys,zs,Ts,ps_quadrupole,'Quadrupole');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: I'm using this snippet as an exact solution to build the true
%       acoustic sources by physical components.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
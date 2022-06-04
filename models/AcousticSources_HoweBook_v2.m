%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Monopole, Dipole and Quadrupoles Generator v2.0
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

% Monopolar Functions
% Ref [1], Eq. 1.7.2
Q = @(t) 1.0*sin(2*pi*t);
psi_monopole = @(r,t) -Q(t-r/c0)./(4*pi*r);

% Dipolar Functions (the most general formulation)
% Ref [1], Eq. 1.7.4, using direct integration of convolution
Ax=1.0;
Ay=1.0;
Az=0.0;
F_1=@(t) Ax*sin(2*pi*t);  dF_1=@(t) Ax*cos(2*pi*t);
F_2=@(t) Ay*sin(2*pi*t);  dF_2=@(t) Ay*cos(2*pi*t);
F_3=@(t) Az*sin(2*pi*t);  dF_3=@(t) Az*cos(2*pi*t);
psi_dipole = @(x1,x2,x3,r,t) 1./(4*pi*c0*r.^3) .* (...
    -x1.*( r.*dF_1(t-r/c0) + c0*F_1(t-r/c0) ) + ...
    -x2.*( r.*dF_2(t-r/c0) + c0*F_2(t-r/c0) ) + ...
    -x3.*( r.*dF_3(t-r/c0) + c0*F_3(t-r/c0) ) ); 

% Lateral Quadrupole Functions (the most general formulation)
% Ref [1], Eq. 1.7.4, using direct integration of convolution
A11=0.0; A12=1.0; A13=1.0;
A21=1.0; A22=0.0; A23=1.0;
A31=1.0; A32=1.0; A33=0.0;
  T_11=@(t)  A11*sin(2*pi*t);   T_12=@(t)  A12*sin(2*pi*t);   T_13=@(t)  A13*sin(2*pi*t);
  T_21=@(t)  A21*sin(2*pi*t);   T_22=@(t)  A22*sin(2*pi*t);   T_23=@(t)  A23*sin(2*pi*t);
  T_31=@(t)  A31*sin(2*pi*t);   T_32=@(t)  A32*sin(2*pi*t);   T_33=@(t)  A33*sin(2*pi*t);
 dT_11=@(t)  A11*cos(2*pi*t);  dT_12=@(t)  A12*cos(2*pi*t);  dT_13=@(t)  A13*cos(2*pi*t);
 dT_21=@(t)  A21*cos(2*pi*t);  dT_22=@(t)  A22*cos(2*pi*t);  dT_23=@(t)  A23*cos(2*pi*t);
 dT_31=@(t)  A31*cos(2*pi*t);  dT_32=@(t)  A32*cos(2*pi*t);  dT_33=@(t)  A33*cos(2*pi*t);
d2T_11=@(t) -A11*sin(2*pi*t); d2T_12=@(t) -A12*sin(2*pi*t); d2T_13=@(t) -A13*sin(2*pi*t);
d2T_21=@(t) -A21*sin(2*pi*t); d2T_22=@(t) -A22*sin(2*pi*t); d2T_23=@(t) -A23*sin(2*pi*t);
d2T_31=@(t) -A31*sin(2*pi*t); d2T_32=@(t) -A32*sin(2*pi*t); d2T_33=@(t) -A33*sin(2*pi*t);
psi_quadrupole = @(x1,x2,x3,r,t) 1./(4*pi*c0^2*r.^5) .* (...
    -x1.*x1.*( r.*d2T_11(t-r/c0) + 3*c0*(r.^2).*dT_11(t-r/c0) + 3*(c0^2)*T_11(t-r/c0)) + ...
    -x1.*x2.*( r.*d2T_12(t-r/c0) + 3*c0*(r.^2).*dT_12(t-r/c0) + 3*(c0^2)*T_12(t-r/c0)) + ...
    -x1.*x3.*( r.*d2T_13(t-r/c0) + 3*c0*(r.^2).*dT_13(t-r/c0) + 3*(c0^2)*T_13(t-r/c0)) + ...
    -x2.*x1.*( r.*d2T_21(t-r/c0) + 3*c0*(r.^2).*dT_21(t-r/c0) + 3*(c0^2)*T_21(t-r/c0)) + ...
    -x2.*x2.*( r.*d2T_22(t-r/c0) + 3*c0*(r.^2).*dT_22(t-r/c0) + 3*(c0^2)*T_22(t-r/c0)) + ...
    -x2.*x3.*( r.*d2T_23(t-r/c0) + 3*c0*(r.^2).*dT_23(t-r/c0) + 3*(c0^2)*T_23(t-r/c0)) + ...
    -x3.*x1.*( r.*d2T_31(t-r/c0) + 3*c0*(r.^2).*dT_31(t-r/c0) + 3*(c0^2)*T_31(t-r/c0)) + ...
    -x3.*x2.*( r.*d2T_32(t-r/c0) + 3*c0*(r.^2).*dT_32(t-r/c0) + 3*(c0^2)*T_32(t-r/c0)) + ...
    -x3.*x3.*( r.*d2T_33(t-r/c0) + 3*c0*(r.^2).*dT_33(t-r/c0) + 3*(c0^2)*T_33(t-r/c0)) );

% Numerical mesh
xa =-5; xb = 5;
ya =-4; yb = 4;
za =-3; zb = 3;
[x,y,z]=meshgrid(xa:0.1:xb,ya:0.1:yb,za:0.1:zb);

% Auxiliary variables
r=sqrt((x-xo).^2+(y-yo).^2+(z-zo).^2);
theta = atan2( (y-yo),(x-xo) ); % azimutal;
phi = atan2( (z-zo),sqrt((x-xo).^2+(y-yo).^2) ); % elevation;

% Compute acoustic pressure fields for time:
p_monopole = psi_monopole(r,tEnd);
p_dipole = psi_dipole(x,y,z,r,tEnd);
p_quadrupole = psi_quadrupole(x,y,z,r,tEnd);

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
    ps_dipole(i,:)     = psi_dipole(xs,ys,zs,r_s,t(i));
    ps_quadrupole(i,:) = psi_quadrupole(xs,ys,zs,r_s,t(i));
end

%% Visualization 2d
figPath = '../figures/';
F1=figure(1); colormap whitejet; plot_source2d(x,y,z, p_monopole ,xs,ys,zs,Ts, ps_monopole , 'Monopole' ,[-0.2,0.2]);
F2=figure(2); colormap whitejet; plot_source2d(x,y,z,  p_dipole  ,xs,ys,zs,Ts,  ps_dipole  ,  'Dipole'  ,[-0.2,0.2]);
F3=figure(3); colormap whitejet; plot_source2d(x,y,z,p_quadrupole,xs,ys,zs,Ts,ps_quadrupole,'Quadrupole',[-0.2,0.2]);

print(F1,[figPath,'directivity_monopole'],'-dpng');
print(F2,[figPath,'directivity_dipole'],'-dpng');
print(F2,[figPath,'directivity_quadrupole'],'-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: These are the most general formulations of the acoustic sources I
%       have obtained. Althought the quadrupole is expensive to compute
%       (and very challenging to obtain), it is very flexible. MD 05/2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
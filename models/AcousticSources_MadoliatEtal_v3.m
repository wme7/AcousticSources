%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Monopole, Dipole and Quadrupoles Generator v3.0
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

% Wave number
k = 9; % defines: k = 2*pi/lambda

% Source location
xo = 0;
yo = 0;
zo = 0;

% Numerical mesh
xa =-5; xb = 5;
ya =-4; yb = 4;
za =-3; zb = 3;
[x,y,z] = meshgrid(xa:0.1:xb,ya:0.1:yb,za:0.1:zb);  

% Dependent variables
r = sqrt((x-xo).^2+(y-yo).^2+(z-zo).^2); 
theta = atan2( (y-yo),(x-xo) ); % azimutal
psi = atan2( (z-zo),sqrt((x-xo).^2+(y-yo).^2) ); % elevation
beta = asin(sin(theta).*cos(psi));

% Monopole
p_monopole  = exp(1i*k*r)./(4*pi*r);
p_dipole    = exp(1i*k*r)./(4*pi*r) .* ((1i*k*r-1)./r) .* cos(theta);
p_quadrupole_long = exp(1i*k*r)./(4*pi.*r.^3) .* ( (3-3*1i*k*r-k*k.*r.^2) .* cos(theta).^2 + 1 - 1i*k*r);
p_quadrupole_latr = exp(1i*k*r)./(4*pi.*r.^3) .* (3-3*1i*k*r-k*k.*r.^2) .* cos(theta) .* cos(beta);

%% Visualization 2d
figure(1); colormap whitejet; plot_source2d(x,y,z, real(p_monopole) , 'Monopole' );
figure(2); colormap whitejet; plot_source2d(x,y,z,  real(p_dipole)  ,  ' Dipole'  );
figure(3); colormap whitejet; plot_source2d(x,y,z,real(p_quadrupole_long),' Quadrupole');
figure(4); colormap whitejet; plot_source2d(x,y,z,real(p_quadrupole_latr),' Quadrupole');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Sadly the source definitions proposed by Madoliat et al (2016) are
% only valid in frequecy space. Furthermore, the lateral quadrupole is not
% satisfactory. (I think there is a mistake/typo). I'll comeback in the
% future to revise these definitions. MD 05/2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Load 
addpath ./third_party/ % third-party tools
addpath ./functions/ % sensors and ploting functions

% Control Parameters
PLOT3d = false;

% Media Properties
c0 = 1.0;

% Source location
xo = 0.0;
yo = 0.0;
zo = 0.0;

% Monopolar Functions
Q = @(t) 1.0*sin(2*pi*t);
psi_monopole = @(r,t) -Q(t-r/c0)./(4*pi*r);

% Dipolar Functions (the general formulation)
F_1 = @(t) 1.0*sin(2*pi*t);
F_2 = @(t) 0.0*sin(2*pi*t);
F_3 = @(t) 0.0*sin(2*pi*t);
psi_dipole = @(x,y,z,r,t) divergence(x,y,z,...
    F_1(t-r/c0)./(4*pi*r),...
    F_2(t-r/c0)./(4*pi*r),...
    F_3(t-r/c0)./(4*pi*r));

% Dipolar Functions (the x_1-directional dipole)
dFdt=@(t) 1.0*cos(2*pi*t); % becuase F := 1.0*sin(2*pi*t);
psi_dipole_2 = @(r,theta,t) dFdt(t-r/c0) .* (cos(theta))./(4*pi*r*c0);

% Lateral Quadrupole Functions (the general formulation)
T_11=@(t) 0.0*sin(2*pi*t); T_12=@(t) 1.0*sin(2*pi*t); T_13=@(t) 1.0*sin(2*pi*t);
T_21=@(t) 1.0*sin(2*pi*t); T_22=@(t) 0.0*sin(2*pi*t); T_23=@(t) 1.0*sin(2*pi*t);
T_31=@(t) 1.0*sin(2*pi*t); T_32=@(t) 1.0*sin(2*pi*t); T_33=@(t) 0.0*sin(2*pi*t);

% Lateral Quadrupole Functions (the (x_1,x_2)-directional quadrupole)
dT2dt2=@(t) -1.0*sin(2*pi*t); % because T := 1.0*sin(2*pi*t);
psi_quadrupole_2 = @(r,theta,phi,t) dT2dt2(t-r/c0) .* (sin(2*theta).*cos(phi))./(8*pi*r*c0^2); 

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

% Compute acoustic pressure fields for time:
t_out = 1.0;
p_monopole = psi_monopole(r,t_out);
%p_dipole = psi_dipole(x,y,z,r,t_out);
p_dipole = psi_dipole_2(r,theta,t_out);
p_quadrupole = psi_quadrupole_2(r,theta,phi,t_out);

%% Build sensors around source
[xs,ys,zs,Ts] = sensors_circle(xo,yo,zo,3.5,0,0,0,100);

% Auxiliary variables for sensors
r_s=sqrt((xs-xo).^2+(ys-yo).^2+(zs-zo).^2);
theta_s = atan2( (ys-yo),(xs-xo) ); % azimulat;
phi_s = atan2( (zs-zo),sqrt((xs-xo).^2+(ys-yo).^2) ); % elevation;

% Sample the pressure field from time t = [0,5]
t = 0:0.1:10; Nt=numel(t); Ns=numel(xs);
ps_monopole   =zeros(Nt,Ns);
ps_dipole     =zeros(Nt,Ns);
ps_quadrupole =zeros(Nt,Ns);
for i = 1:Nt
    ps_monopole(i,:)   = psi_monopole(r_s,t(i));
    ps_dipole(i,:)     = psi_dipole_2(r_s,theta_s,t(i));
    ps_quadrupole(i,:) = psi_quadrupole_2(r_s,theta_s,phi_s,t(i));
end

%% Visualization 2d
figure(1); colormap whitejet;
    subplot(221)
        sm=surf(x(:,:,zc),y(:,:,zc),p_monopole(:,:,zc)); hold on
        scatter3(xs,ys,zs,'.k'); hold off
        axis equal; title('Monopole'); shading interp;
        xlabel('x_1'); ylabel('x_2'); view(2); colorbar; caxis([-0.4,0.4]);
    subplot(223)
        imagesc(ps_monopole);
        title('Signal enregistré'); xlabel('capteurs'); ylabel('time [ms]');
    subplot(2,2,[2,4])
        polarplot(theta2,rms(ps_monopole)); thetatickformat('degrees')
figure(2); colormap whitejet;
    subplot(221)
        sd=surf(x(:,:,zc),y(:,:,zc),p_dipole(:,:,zc)); hold on
        scatter3(xs,ys,zs,'.k'); hold off
        axis equal; title('Dipole'); shading interp;
        xlabel('x_1'); ylabel('x_2'); view(2); colorbar; caxis([-0.2,0.2]);
    subplot(223)
        imagesc(ps_dipole);
        title('Signal enregistré'); xlabel('capteurs'); ylabel('time [ms]');
    subplot(2,2,[2,4])
        polarplot(theta2,rms(ps_dipole)); thetatickformat('degrees')
figure(3); colormap whitejet;
    subplot(221)
        sq=surf(x(:,:,zc),y(:,:,zc),p_quadrupole(:,:,zc)); hold on
        scatter3(xs,ys,zs,'.k'); hold off
        axis equal; title('Quadrupole'); shading interp;
        xlabel('x_1'); ylabel('x_2'); view(2); colorbar; caxis([-0.1,0.1]);
    subplot(223)
        imagesc(ps_quadrupole);
        title('Signal enregistré'); xlabel('capteurs'); ylabel('time [ms]');
    subplot(2,2,[2,4])
        polarplot(theta2,rms(ps_quadrupole)); thetatickformat('degrees')

%% Visualization 3d
if PLOT3d
figure(4); colormap jet;
    isosurface(x,y,z,p_monopole);
    axis equal; title('Monopole'); alpha(0.5); 
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); view(3);
figure(5); colormap jet;
    isosurface(x,y,z,p_dipole);
    axis equal; title('Dipole'); alpha(0.5); 
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); view(3);
figure(6); colormap jet;
    isosurface(x,y,z,p_quadrupole);
    axis equal; title('Quadrupole'); alpha(0.5); 
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); view(3);
end

%% Extra functions
function [x_i,y_i,z_i,theta_i] = sensors_circle(xo,yo,zo,radius,yaw,pitch,roll,nSensors)
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
% Integrate pressure field
clear; close all;

% Load 
addpath ../third_party/ % third-party tools
addpath ../functions/ % sensors and ploting functions

% Define model parameters
L = 1; k = 6;

% Observers coordinates
r =@(x1,x2) sqrt(x1.^2+x2.^2);
T =@(x1,x2) atan2(x2,x1);
z =@(x3)    x3;

% Source coordinates
ro =@(y1,y2) sqrt(y1.^2+y2.^2);
To =@(y1,y2) atan2(y2,y1);
zo =@(y3)    y3;

% Wave packet function
Txx =@(x) exp(-x.^2/L^2) .* exp(1i*k*x);
dxx_Txx =@(x) ((-2*x/L^2 + 1i*k).^2 - 2/L^2) .* Txx(x);

% Define free Stream functions
Go =@(R) exp(-1i*k*R)./R;
%dxx_Go = d2GodR2.*dRdx.^2 + dGodR.*d2Rdx2;

% Define distances R
R =@(r,T,z,ro,To,zo) sqrt(r.^2+ro.^2-2*r.*ro.*cos(T-To)+(z-zo).^2);

% Integrate pressure field
%Kernel1 =@(x1,x2,x3,y1,y2,y3) Txx(y1).*dxx_Go(...
%    R(r(x1,x2),T(x1,x2),z(x3),ro(y1,y2),To(y1,y2),zo(y3)));
Kernel2 =@(x1,x2,x3,y1,y2,y3) dxx_Txx(y1).*Go(...
    R(r(x1,x2),T(x1,x2),z(x3),ro(y1,y2),To(y1,y2),zo(y3)));

% Define mesh points (observer points)
[x1,x2,x3,Ts] = build_circularArrayOfSensors(0,0,0,35*L,60);
theta = Ts*(180/pi);

% Compute pressure fields
pressure1 = zeros(size(x1));
pressure2 = zeros(size(x1));

parfor i=1:numel(x1)
    % Define function to be integrated
    %f1 =@(y1,y2,y3) Kernel1(x1(i),x2(i),x3(i),y1,y2,y3);
    f2 =@(y1,y2,y3) Kernel2(x1(i),x2(i),x3(i),y1,y2,y3);

    % The pressure field
    %pressure1(i) = 1/(4*pi)*integral3(f1,-12,12,-12,12,-12,12);
    pressure2(i) = 1/(4*pi)*integral3(f2,-12,12,-12,12,-12,12);
end
save('_sensors.mat','theta','pressure2');

% Visualize sensors
%plot(theta,real(pressure1),'-'); hold on
plot(theta,real(pressure2),'o'); hold off
xlabel('Polar'); ylabel('SPL');

% Visualize fields
%slice(x1,x2,x3,real(pressure1),[],0,0); shading interp; colormap whitejet;
%slice(x1,x2,x3,real(pressure2),[],0,0); shading interp; colormap whitejet;
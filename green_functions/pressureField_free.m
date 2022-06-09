% Integrate pressure field
clear; close all;

% Load 
addpath ../third_party/ % third-party tools
addpath ../functions/ % sensors and ploting functions

% Define model parameters
L = 1;   % Some charac length
U = 0.6; % also the Mach number
St= 0.2; % St = f*L/U
f = St*U/L; % frequency
k = 2*pi*f; % wave number
c = 1.0;  % speed of sound
omega = c*k; % angular frequency

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

% Define distances R
R =@(r,T,z,ro,To,zo) sqrt(r.^2+ro.^2-2*r.*ro.*cos(T-To)+(z-zo).^2);

% dxx_Go = d2GodR2.*dRdx.^2 + dGodR.*d2Rdx2;
dGodR =@(R) -(1./R+ 1i*k).*Go(R);
d2GodR2 =@(R) -(k+2./(R.^2)+ 2i*k./R).*Go(R);
dRdro =@(R,r,T,ro,To) (ro-r.*cos(T-To))./R;
dRdTo =@(R,r,T,ro,To) -ro.*r.*sin(T-To)./R;
dRdx =@(R,r,T,ro,To) cos(To).*dRdro(R,r,T,ro,To) - ...
    (sin(To)./ro).*dRdTo(R,r,T,ro,To);
d2Rdro2 =@(R,r,T,ro,To) 1./R - (ro-r.*cos(T-To))./(2*R.^3);
d2RdTo2 =@(R,r,T,ro,To) r.*ro.*cos(T-To)./R - ((r.*ro.*sin(T-To)).^2)./R.^3;
d2Rdx2 =@(R,r,T,ro,To) (cos(To).^2) .* d2Rdro2(R,r,T,ro,To) - ...
    2*sin(To).* cos(To).* dRdro(R,r,T,ro,To).* dRdTo(R,r,T,ro,To)./ro + ...
    (sin(To).^2).* d2RdTo2(R,r,T,ro,To)./ro;
dxx_Go =@(R,r,T,ro,To) d2GodR2(R).*dRdx(R,r,T,ro,To).^2 + ...
         dGodR(R).*d2Rdx2(R,r,T,ro,To);

% Integrate pressure field
Kernel1 =@(x1,x2,x3,y1,y2,y3) Txx(y1).*dxx_Go(...
	R(r(x1,x2),T(x1,x2),z(x3),ro(y1,y2),To(y1,y2),zo(y3)),...
    r(x1,x2),T(x1,x2),ro(x1,x2),To(x1,x2));
Kernel2 =@(x1,x2,x3,y1,y2,y3) dxx_Txx(y1).*Go(...
    R(r(x1,x2),T(x1,x2),z(x3),ro(y1,y2),To(y1,y2),zo(y3)));

% Define mesh points (observer points)
nSensors = 5;
[x1,x2,x3,Ts] = build_circularArrayOfSensors(0,0,0,35*L,nSensors);
theta = Ts*(180/pi);

% Compute pressure fields
pressure1 = zeros(size(x1));
pressure2 = zeros(size(x1));

%% Numerical Integration
dy=0.1;
[y1,y2,y3] = meshgrid(-7:dy:7);

% Prepare figures
array1 = zeros(size(y1));
figure(1); colormap whitejet;
figure(2); colormap whitejet;

for i=1:nSensors
    % Define function to be integrated
    %f1 =@(y1,y2,y3) Kernel1(x1(i),x2(i),x3(i),y1,y2,y3);
    %f2 =@(y1,y2,y3) Kernel2(x1(i),x2(i),x3(i),y1,y2,y3);
    array1 = Kernel1(x1(i),x2(i),x3(i),y1,y2,y3);
    array2 = Kernel2(x1(i),x2(i),x3(i),y1,y2,y3);
    
    % Visualize fields
    F1=figure(1); slice(y1,y2,y3,real(array1),[],0,0); shading interp;
    F2=figure(2); slice(y1,y2,y3,real(array2),[],0,0); shading interp;
    drawnow;

    % Integrate the pressure field
    %pressure1(i) = 1/(4*pi)*integral3(f1,-12,12,-12,12,-12,12);
    %pressure2(i) = 1/(4*pi)*integral3(f2,-12,12,-12,12,-12,12);
    pressure1(i) = 1/(4*pi)*trapz3d(array1,dy,dy,dy);
    pressure2(i) = 1/(4*pi)*trapz3d(array2,dy,dy,dy);
end
pressure1(i+1)=pressure1(1);
pressure2(i+1)=pressure2(1);
save('_sensors_FREE_.mat','theta','pressure1','pressure2');

%% Visualize sensors
F3=figure(3);
figPath = '../figures/';
plot(theta,real(pressure1),'-x'); hold on
plot(theta,real(pressure2),'-o'); hold off
xticks(-180:30:180);
xlabel('$\theta_{Polar}$','Interpreter','latex','Fontsize',20);
ylabel('SPL (arbitrary dB)','Interpreter','latex','Fontsize',20); 
title('Free; M=0.6; St=0.2','Interpreter','latex','Fontsize',20);
legend({'$\partial^2 T_{xx}/\partial x^2$','$\partial^2 G_o/\partial x^2$'},...
    'location','best','Interpreter','latex','Fontsize',20);

print(F1,[figPath,'pressureField_FREE_method1'],'-dpng');
print(F2,[figPath,'pressureField_FREE_method2'],'-dpng');
print(F3,[figPath,'pressureField_FREE_compare'],'-dpng');
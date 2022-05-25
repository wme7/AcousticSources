% Sources V2

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

% Wave number
k = 9; % defines: k = 2*pi/lambda

% Monopole
p_monopole  = exp(1i*k*r)./(4*pi*r);
p_dipole    = exp(1i*k*r)./(4*pi*r) .* ((1i*k*r-1)./r) .* cos(theta);
p_quad_long = exp(1i*k*r)./(4*pi.*r.^3) .* ( (3-3*1i*k*r-k*k.*r.^2) .* cos(theta).^2 + 1 - 1i*k*r);
p_quad_latr = exp(1i*k*r)./(4*pi.*r.^3) .* (3-3*1i*k*r-k*k.*r.^2) .* sin(2*theta) .* cos(beta);

%% Visualization
figure(1); isosurface(x,y,z,real(p_monopole)); title('monopole'); alpha(0.5);
figure(2); isosurface(x,y,z,real(p_dipole)); title('dipole'); alpha(0.5);
figure(3); isosurface(x,y,z,real(p_quad_long)); title('quadrupole - longitudinal'); alpha(0.5);
figure(4); isosurface(x,y,z,real(p_quad_latr)); title('quadrupole - lateral'); alpha(0.5);
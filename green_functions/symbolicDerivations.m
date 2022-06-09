% Produce symbolic derivations
clear, clc;

% Properties
syms k

% Cylindrical coords
syms r ro z zo T To

% Cartesian coords
syms x1 x2 x3

% Define distances R and R':
R = sqrt(r.^2+ro.^2-2*r.*ro.*cos(T-To)+(z-zo).^2);
R_= sqrt(r.^2+ro.^2-2*r.*ro.*cos(T+To)+(z-zo).^2);

% Derivatives of R:
dRdro = simplify(expand(diff(R,ro)));
dRdTo = simplify(expand(diff(R,To)));
d2Rdro2 = simplify(expand(diff(R,ro,2)));
d2RdTo2 = simplify(expand(diff(R,To,2)));

% Define B
B = (r+ro).^2+(z-zo).^2;

% Derivatives of B:
dBdro = simplify(expand(diff(B,ro)));
d2Bdro2 = simplify(expand(diff(B,ro,2)));

%clear; syms r ro z zo T To R R_ B k
% Define distances uR and uR':
uR = 2*sqrt(k*r*ro./( B+R )).*cos(0.5*(T-To));
uR_= 2*sqrt(k*r*ro./( B+R_)).*cos(0.5*(T+To));

% Derivatives of uR:
%duRdR = simplify(expand(diff(uR,R)));
%duRdB = simplify(expand(diff(uR,B)));
duRdro = simplify(expand(diff(uR,ro)));
duRdTo = simplify(expand(diff(uR,To)));
d2uRdro2 = simplify(expand(diff(uR,ro,2)));
d2uRdTo2 = simplify(expand(diff(uR,To,2)));

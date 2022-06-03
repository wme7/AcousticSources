% ===============================================================
%     Code used only to generate and save the integral tables
% ===============================================================
function generateTables

% Generate the integral tables, more accurate than Abramowitz &
% Stegun provide, since they give only 7 digits.
FresnelCObj = @(t) cos(pi*t.^2/2);
FresnelSObj = @(t) sin(pi*t.^2/2);

p = 1.75;
T0 = linspace(1,7.5.^p,501).' .^(1/p);
dt = T0(2) - T0(1);
T0 = [linspace(0,1 - dt,ceil(1./dt))';T0];
plot(diff(T0))

n = length(T0);
FC75 = zeros(n,1);
FS75 = zeros(n,1);

h = waitbar(0,'Computing Fresnel integrals');
for i = 2:n
  waitbar(i/n,h)
  FC75(i) = quadgk(FresnelCObj,0,T0(i),'abstol',1.e-16,'reltol',100*eps('double'));
  FS75(i) = quadgk(FresnelSObj,0,T0(i),'abstol',1.e-16,'reltol',100*eps('double'));
end
delete(h)

% Turn them into splines, then save the splines. These splines are
% first built in a Hermite form, since I can supply the 1st and second
% derivatives of the function. Then I turn them into a pp form, for use
% in fresnelC and fresnelS.
FCspl = hermite2slm([T0,FC75,FresnelCObj(T0), -pi*T0.*sin(pi*T0.^2/2), ...
  -pi*(sin(pi*T0.^2/2) + pi*T0.^2 .*cos(pi*T0.^2/2))]);
FCspl = slm2pp(FCspl);

FSspl = hermite2slm([T0,FS75,FresnelSObj(T0),pi*T0.*cos(pi*T0.^2/2), ...
  pi*(cos(pi*T0.^2/2) - pi*T0.^2 .*sin(pi*T0.^2/2))]);
FSspl = slm2pp(FSspl);

save _Fresnel_data_ FCspl FSspl


% test the result
clear functions

n = 1000;
T = sort(rand(n,1)*10);
FCquad = zeros(n,1);
FSquad = zeros(n,1);
for i = 1:n
  FCquad(i) = quadgk(FresnelCObj,0,T(i),'abstol',1.e-16);
  FSquad(i) = quadgk(FresnelSObj,0,T(i),'abstol',1.e-16);
end
FCpred = fresnelC(T,0);
FSpred = fresnelS(T,0);

subplot(1,2,1)
plot(T,FCquad - FCpred,'.')
grid on
subplot(1,2,2)
plot(T,FSquad - FSpred,'.')
grid on

end
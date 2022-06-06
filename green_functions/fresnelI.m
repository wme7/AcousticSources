function [IR] = fresnelI(u)
% fresnelI: computes the complex integral of exp(-i * t^2) from -inf to u.
%
%             Coded by Manuel A. Diaz, Univ-Poitiers, 2022.
%
% Assumptions: 
% (1) we use fresnelC and fresnelS to approximate the complex exponential
%     as:  exp(-i*u^2) = cos(u^2) - i*sin(u^2),
% (2) We know that the fresnel intregral from 0 to Inf is equal to 
%     sqrt(pi/2)*(1-i)/2. The negative argument is assumed for -Inf to 0.
%     Ref[1] https://en.wikipedia.org/wiki/Fresnel_integral.
%
IR = sqrt(0.5*pi) * (fresnelC(u,1) - 1i*fresnelS(u,1) + 0.5 - 0.5i);
end


function output = trapz3d(mat,dx,dy,dz)
%
output = trapz(dx,trapz(dy,trapz(dz,mat,3),2),1);
%
end
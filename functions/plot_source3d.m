function plot_source3d(x,y,z,field)
    colormap jet;
    isosurface(x,y,z,field);
    axis equal; title('Monopole'); alpha(0.5); 
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); view(3);
end
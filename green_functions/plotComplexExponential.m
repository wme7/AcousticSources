% Visualize the complex field produced by f(x) = exp(-1i*x);

% Numerical grid
t = 0:0.05:2*pi;

% Wave number
k = 4.0;

% Evaluate the term
complex_field = exp(-1i*k*t);

% Visualize (Real-Imag) vs. time
plot3(t, real(complex_field), imag(complex_field)); grid on;
xlabel('time','interpreter','latex','Fontsize',20);
ylabel('$\Re(e^{-i\,x})$','interpreter','latex','Fontsize',20);
zlabel('$\Im(e^{-i\,x})$','interpreter','latex','Fontsize',20);
xticks(0:pi/4:2*pi);
% Visualize the complex field produced by f(x) = exp(-1i*x);

% Numerical grid
x = 0:0.1:10;

% Evaluate the term
complex_field = exp(1i*x);

% Visualize
plot3(x, real(complex_field), imag(complex_field)); grid on;
ylabel('\Re(e^{1i*x})');
zlabel('\Im(e^{1i*x})');
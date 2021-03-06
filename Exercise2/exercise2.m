% Exercise 1: Statistics and analyzing estimators

deltaT = [5; 15; 35; 60];
L = [6.55; 9.63; 17.24; 29.64];

% (a) Plot the data
hold on;
plot(deltaT, L,"x");    

% (b) Experimental values of L_0 and A:
% WARNING: do not rename variables this will break the tests!
% Computation on Paper
b = [sum(L); sum(deltaT .* L)];
X = [sum(deltaT), size(L, 1) ; sum(deltaT .* deltaT), sum(deltaT)];

res = inv(X) * b;

A = res(1);
L_0 = res(2);

p = polyfit(deltaT, L, 1);
disp(p);

x1 = linspace(-10, 90);
y1 = polyval(p,x1);
plot(x1, y1,"x");

% (c) Fit third order polynomial to data
% WARNING: do not rename variables this will break the tests!
poly_coeffs = polyfit(deltaT, L, 3);
y1 = polyval(poly_coeffs,x1);
plot(x1, y1, "x");

legend('Points', '1 Degree', '3 Degree');
xlabel("Delta T");
ylabel("L");

% (d) Validate fits with additional measurement
deltaT_val = 70;
L_val = 32.89;

plot(deltaT_val, L_val,"x");

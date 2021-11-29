function results = wlsExercise(I, U)

N_e = size(I, 1);    % number of students/experiments
N_m = size(I, 2);    % number of measurements per experiment
disp(N_m);
disp(N_e);

%% Plot data 

figure(1); hold on;

for k=1:N_e
    plot(I(k, :), U(k, :), 'x');
end

xlabel('I'); ylabel('U');


%% WLS/LLS for first student
i1 = I(1, :)';
Phi = [ones(size(i1, 1), 1) i1];

theta_LLS_1 = inv(Phi' * Phi) * Phi' * U(1, :)';
U_LLS_1 = Phi * theta_LLS_1;

c = 1;
var = c * (1:1:N_m);
disp(var);
W = diag(1 ./ var);
disp(size(W));
disp(size(Phi));

theta_WLS_1 = inv(Phi' * W * Phi) * Phi' * W * U(1, :)';
U_WLS_1 = Phi * theta_WLS_1;

figure(2); hold on;
plot(I(1, :), U(1, :), 'x');
plot(I(1, :), U_LLS_1, '-');
plot(I(1,:), U_WLS_1);

legend('data', 'LLS', 'WLS', 'Location', 'northwest');

%% LLS/WLS for all students

thetas_LLS = zeros(N_e, 2);
thetas_WLS = zeros(N_e, 2);

for i=1:N_e
    i1 = I(i, :)';
    Phi = [ones(size(i1, 1), 1) i1];

    thetas_LLS(i, :) = inv(Phi' * Phi) * Phi' * U(i, :)';
    thetas_WLS(i, :) = inv(Phi' * W * Phi) * Phi' * W * U(i, :)';
end
    
%% Estimate covariance of theta

theta_mean_LLS = (1 / N_e) * sum(thetas_LLS);
sigma_LLS = (1 / (N_e - 1)) * (thetas_LLS - theta_mean_LLS)' * (thetas_LLS - theta_mean_LLS);

theta_mean_WLS = (1 / N_e) * sum(thetas_WLS);
sigma_WLS = (1 / (N_e - 1)) * (thetas_WLS - theta_mean_WLS)' * (thetas_WLS - theta_mean_WLS);

% Plot confidence ellipsoids

% Compute the eigenvalues and eigenvectors for estimators
[V_LLS, D_LLS] = eig(sigma_LLS); %eig_values and eig_vectors
[V_WLS, D_WLS] = eig(sigma_WLS); %eig_values and eig_vectors

% Generate coordinates of 50 points on a unit circle
xy = [cos(linspace(0, 2*pi, 50)); sin(linspace(0, 2*pi, 50))];

% Generate the points of the confidence ellipse
xy_ellipse1 = theta_mean_LLS'*ones(1,50) + V_LLS*sqrt(D_LLS)*xy;
xy_ellipse2 = theta_mean_WLS'*ones(1,50) + V_WLS*sqrt(D_WLS)*xy;

figure(3); hold on;
plot(thetas_LLS(:, 1), thetas_LLS(:, 2), 'rx');
plot(thetas_WLS(:, 1), thetas_WLS(:, 2), 'bx');
plot(xy_ellipse1(1, :), xy_ellipse1(2, :), 'r-');
plot(xy_ellipse2(1, :), xy_ellipse2(2, :), 'b-');

plot(theta_mean_LLS(1), theta_mean_LLS(2), 'r.', 'MarkerSize', 20);
plot(theta_mean_WLS(1), theta_mean_WLS(2), 'b.', 'MarkerSize', 20);

legend('LLS', 'WLS', 'LLS', 'WLS', 'LLS', 'WLS');

%% return all relevant results as a struct

results = struct('theta_LLS_1', theta_LLS_1, ...
                 'theta_WLS_1', theta_WLS_1, ...
                 'U_LLS_1', U_LLS_1, ...
                 'U_WLS_1', U_WLS_1, ...
                 'thetas_LLS', thetas_LLS, ...
                 'thetas_WLS', thetas_WLS, ...
                 'theta_mean_LLS', theta_mean_LLS, ...
                 'theta_mean_WLS', theta_mean_WLS, ...
                 'sigma_LLS', sigma_LLS, ...
                 'sigma_WLS', sigma_WLS);

end

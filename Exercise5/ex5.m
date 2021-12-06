function result = ex5(X, Y)

%% Setup
alphas = [0, 10e-6, 10e-5, 1];

N_alpha = length(alphas);

order = 7;

N_theta = order+1;   % number of parameters

N_e = size(X, 1);    % number of students/experiments
N_m = size(X, 2);    % number of measurements per experiment

%% x and Phi used for plots

xs = linspace(0,1, 200);
PhiTest = ones(200, order+1);

for o=1:order
    PhiTest(:, o+1) = xs'.^o;
end

%% Single experiment

x1 = X(1, :)';
y1 = Y(1, :)';
Phi = [ones(N_m, 1) x1 x1.^2 x1.^3 x1.^4 x1.^5 x1.^6 x1.^7];

thetas_single = zeros(N_alpha, N_theta);
thetas_single_norm = zeros(N_alpha, 1);
R_squared = zeros(N_alpha, 1);

figure(1); hold on; grid on;

plot(X(1, :), Y(1, :), 'x');
xlim([0,1]); ylim([-7,7]);
xlabel('X'); ylabel('Y');

for a=1:N_alpha
    alphaI = alphas(a) * eye(N_theta);
    
    thetas_single(a, :) = inv(Phi' * Phi + alphaI) * Phi' * y1;
    theta = thetas_single(a, :)';
    
    thetas_single_norm(a) = norm(theta);
%     R_squared(a) = (norm(Phi * theta) / norm(y1)) ^ 2;
    R_squared(a) = 1 - (norm(y1 - Phi * theta) ^ 2 / norm(y1) ^ 2);
    
    plot(xs, PhiTest*theta, '-');
end

disp(R_squared);

legend('data', 'no regularization', 'small regularization', 'strong regularization', ...
       'very strong regularization', 'Location', 'southwest');

%% Multiple Experiments

thetas = zeros(N_alpha, N_e, N_theta);
thetas_mean = zeros(N_alpha, N_theta);
    
figure(2);
titles = {'no regularization', 'small regularization', ...
          'strong regularization', 'very strong'};

for a=1:N_alpha
    
    alphaI = alphas(a) * eye(N_theta);
        
    subplot(N_alpha, 2, 2*a - 1); hold on; grid on;
    title(titles{a});
    ylim([-6, 6]); xlim([0, 1]);
    
    for k=1:N_e
        xk = X(k, :)';
        yk = Y(k, :)';

        Phi = [ones(N_m, 1) xk xk.^2 xk.^3 xk.^4 xk.^5 xk.^6 xk.^7];

        theta = inv(Phi' * Phi + alphaI) * Phi' * yk;
        thetas(a, k, :) = theta';
        plot(xs, PhiTest*theta, '-');  
    end
    
    subplot(N_alpha, 2, 2*a); hold on; grid on;
    title(titles{a});
    ylim([-6, 6]); xlim([0, 1]);
    
    theta_mean = squeeze(mean(thetas(a, :, :), 2));
    disp(size(theta_mean));

    thetas_mean(a, :) = theta_mean;
    
    plot(xs, PhiTest*theta_mean, '-');
end

result = struct('thetas_single', thetas_single, ...
                'thetas_single_norm', thetas_single_norm, ...
                'R_squared', R_squared, ...
                'thetas', thetas, ...
                'thetas_mean', thetas_mean);
end

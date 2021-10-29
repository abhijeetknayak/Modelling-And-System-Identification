N = 1000; % number of samples
M = 200;  % number of experiments

% load exercise data
load('exercise1_dataset.mat'); % provides: i_k (currents) and v_k (voltages)

%% (a) 
% Compute the following three quantities 
% WARNING: do not rename variables this will break the tests!

R_SA_single = zeros(N, 1);
R_LS_single = zeros(N, 1);
R_EV_single = zeros(N, 1);

for idx = 1:N 
    R_SA_single(idx, 1) = mean(v_k(1:idx, 1) ./ i_k(1:idx, 1));
    R_LS_single(idx, 1) = mean(v_k(1:idx, 1) .* i_k(1:idx, 1)) / mean(i_k(1:idx, 1).* i_k(1:idx, 1));
    R_EV_single(idx, 1) = mean(v_k(1:idx, 1)) / mean(i_k(1:idx, 1));
end

% plot the computed quantities in a single plot
figure(1); % create your plots in this figure
hold on;
title("Resistance Single Student");
xlabel("Measurement");
ylabel("Resistance");
plot(1:N, R_SA_single, Marker="*", Color='b');
plot(1:N, R_LS_single, Marker=".", Color='g');
plot(1:N, R_EV_single, Marker="+", Color='r');

legend("SA", "LS", "EV");

%% (b)
% WARNING: do not rename variables this will break the tests!

R_SA = zeros(N, M);
R_LS = zeros(N, M);
R_EV = zeros(N, M);

for student=1:M
    for idx = 1:N 
        R_SA(idx, student) = mean(v_k(1:idx, student) ./ i_k(1:idx, student));
        R_LS(idx, student) = mean(v_k(1:idx, student) .* i_k(1:idx, student)) / mean(i_k(1:idx, student).* i_k(1:idx, student));
        R_EV(idx, student) = mean(v_k(1:idx, student)) / mean(i_k(1:idx, student));
    end
end

% plots
figure(2);
hold on;
title("Simple Approach");
xlabel("Measurements");
ylabel("Resistance");
for idx=1:M
    plot(1:N, R_SA(:, idx), Marker="*", Color='b',LineWidth=0.025);
end

figure(3);
hold on;
title("Least Squares");
xlabel("Measurements");
ylabel("Resistance");
for idx=1:M
    plot(1:N, R_LS(:, idx), Marker=".", Color='g',LineWidth=0.025);
end

figure(4);
hold on;
title("Error in Variables");
xlabel("Measurements");
ylabel("Resistance");
for idx=1:M
    plot(1:N, R_EV(:, idx), Marker="+", Color='r',LineWidth=0.025);
end


%% (c)
% WARNING: do not rename variables this will break the tests!
R_SA_mean = mean(R_SA, 2);
R_LS_mean = mean(R_LS, 2);
R_EV_mean = mean(R_EV, 2);

% plot
figure(5);
hold on
title("Resistance Mean");
xlabel("Measurement");
ylabel("Resistance");

plot(1:N, R_SA_mean, Marker="*", Color='b',LineWidth=0.025);
plot(1:N, R_LS_mean, Marker=".", Color='g',LineWidth=0.025);
plot(1:N, R_EV_mean, Marker="+", Color='r',LineWidth=0.025);

legend("SA", "LS", "EV")

%% (d)
% WARNING: do not rename variables this will break the tests!
R_SA_Nmax = R_SA(N, :);
R_LS_Nmax = R_LS(N, :);
R_EV_Nmax = R_EV(N, :);

%plots
figure(6);
histogram(R_SA_Nmax, 100);
title("Simple Approach");
xlabel("Resistance");
ylabel("Occurences");

figure(7);
histogram(R_LS_Nmax, 20);
title("Least Squares");
xlabel("Resistance");
ylabel("Occurences");

figure(8);
histogram(R_EV_Nmax, 20);
title("Error in Variables");
xlabel("Resistance");
ylabel("Occurences");
%% WARNING do not rename any variable, this will break the tests!

%% load data
data = load('exercise3_dataset.mat');
i1 = data.i1;
u1 = data.u1;
i2 = data.i2;
u2 = data.u2;
%% Plot data
% Use subplot to plot both data sets and their linear fits (later in b)) in two separate plots, but in the same plotting window
figure(1)
plot1 = subplot(2,1,1);
% code for first plot
plot(i1, u1);
title("I vs V");
xlabel("I");
ylabel("V");

% code for second plot
hold on
plot2 = subplot(2,1,2);
plot(i2, u2);
title("I vs V");
xlabel("I");
ylabel("V");

%% Compute Linear least squares
Phi1 = [ones(size(i1, 1), 1) i1];
Phi2 = [ones(size(i2, 1), 1) i2];

theta1 = inv(Phi1' * Phi1) * Phi1' * u1;
theta2 = inv(Phi2' * Phi2) * Phi2' * u2;


% % Assign the estimated values
R_est1 = theta1(2);
R_est2 = theta2(2);
E_est1 = theta1(1);
E_est2 = theta2(1);

%Plot linear fit into the arlier created plot
subplot(plot1) %using earlier saved handle to plot in the first subplot
hold on;
% i = linspace(0, 0.01, 1000);
u_pred1 = E_est1 + R_est1 .* i1;
plot(i1, u_pred1);

subplot(plot2); %using earlier saved handle
hold on;
% i = linspace(-0.01, 0.01, 2000);
u_pred2 = E_est2 + R_est2 .* i2;
plot(i2, u_pred2);

% %% Calculate Residuals 
r1 = (u1 - u_pred1);
r2 = (u2 - u_pred2);
% Plot histogram of them
figure(2)
subplot(1,2,1) % use subplot once again
hold on
% code for first histogram
histogram(r1, 100);
xlabel("Value");
ylabel("Count");
title("Histogram for u1");

subplot(1,2,2)
hold on;
% code for second histogram
histogram(r2, 100);
xlabel("Value");
ylabel("Count");
title("Histogram for u2");

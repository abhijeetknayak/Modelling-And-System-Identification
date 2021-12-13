%load data
data = load('exercise6_data.mat');
% wind speed data
v_w = data.v_w;
% Plot histogram
figure(1);
h1 = histogram(v_w);
h1.Normalization = 'probability';
h1.BinWidth = 1.0;
title('Histogram of the wind speeds at Feldberg');
ylabel('frequency of occurence');
xlabel('v_w [m/s]');
    
%% Estimation of lambda and K of Weibull distribution
% Generate wind samples
N = length(v_w);
% give lambda and k initial values for fmincon
lambda_0 = 2;
k_0 = 10;
% define the initial values in a vector for fmincon    
theta0 = [k_0, lambda_0];
b=log(v_w);

% implement the function fun corresponding to the objective function derived in b)
fun = @(x) N*(x(1)*log(x(2))-log(x(1)))-(x(1)-1)*sum(b)+(x(2)^(-x(1)))*sum(v_w.^x(1));

% options for fmincon
% options = optimoptions('fmincon','Display','iter');
% call fmincon to solve the problem
A=[];
b=[];
Aeq=[];
beq=[];
lb= [0,0];
ub=[inf,inf];
[theta,y] = fmincon(fun,theta0,A,b,Aeq,beq,lb,ub);

% define a x-Axis vector to plot the fit
v = linspace(0, 30, 1000);
% implement the probability density function given p(v|k,lambda) in order to plot the fit
pDist = @(v,k,l)(k/l).*((v/l).^(k-1)).*exp(-(v/l).^k);
% %% Plot fitted pfd into histogram
figure(2);
h1 = histogram(v_w);
h1.Normalization = 'probability';
h1.BinWidth = 1.0;
hold on;
plot(v, pDist(v, theta(1), theta(2)), 'r', 'lineWidth', 3);
legend('Histogram of v_w', 'fitted Weibull distribution');
title('PDF of the wind speed');
ylabel('p(v_w)');
xlabel('v_w [m/s]')

%extract the optimal parameters from the vector theta
kStar = theta(1);
lambdaStar = theta(2);
    
%% Expected power
% implement the weibull distribution function with the optimal parameter
Weibull = @(v) (kStar/lambdaStar).*((v/lambdaStar).^(kStar-1)).*exp(-(v/lambdaStar).^kStar)

% declare a sample vector with wind speeds
WindSpeed = linspace(0,26, 27);
% compute the probabilities using the Weibull distribution
windSpeedProbability = Weibull(WindSpeed);
% define a vector containing the power given in the table
power = [0; 0; 3; 25; 82; 174; 321; 532; 815; 1180; 1580; 1900; 2200; 2480; 2700; 2850; 2950; 3020; 3020; 3020; 3020; 3020; 3020; 3020; 3020; 3020; 0];

expectedPower = 0;
% Compute the expected power using the trapezoid rule  
for i=1:26
   expectedPower = expectedPower+(power(i)+power(i+1))/2
end

fprintf('k* = %0.3f | lambda* = %0.3f | E(Power)* = %0.3f \n', kStar, lambdaStar, expectedPower);
               
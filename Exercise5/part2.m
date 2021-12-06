%% LLS Problem Definition
N = 9;
y = (1:N)';
phi = [ones(N,1)/sqrt(N) -ones(N,1)/sqrt(N)];

%% Your code goes here. Warning: Do not rename the variables!

%You can comment this out if you just want to see the minimization problem.
alpha = 0.1;
theta_reg = inv(phi' * phi + alpha * eye(2)) * phi' * y;
theta_MPPI = pinv(phi' * phi) * phi' * y;

%% Visualization
t1 = linspace(5,10,50);
t2 = linspace(-10,-5,50);

Q = phi'*phi; c = -2*y'*phi;
f_visu = @(t1,t2) y'*y + c(1)*t1 +c(2)*t2 + Q(1,1)*t1.^2 +2*Q(1,2).*t1.*t2 + Q(2,2)*t2.^2;
[X1,X2] = meshgrid(t1,t2);
Z = f_visu(X1,X2);

hold on
%plot surface
surf(X1,X2,Z);
%plot solutions
if exist('theta_reg','var')
    plot3(theta_reg(1),theta_reg(2),f_visu(theta_reg(1),theta_reg(2)),'rx');
end
if exist('theta_MPPI','var')
    plot3(theta_MPPI(1),theta_MPPI(2),f_visu(theta_MPPI(1),theta_MPPI(2)),'rx');
end

xlabel('\theta_1');ylabel('\theta_2');zlabel('f(\theta)');view(-22,55);
hold off
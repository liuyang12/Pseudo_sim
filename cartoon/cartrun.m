%CARTRUN run of cartoon model
%% [cartrhs]
N = 5;
I = 0.1;
D = 0.5;
cartfun = @(t,y) cartrhs(t,y,I,D);
tspan = [0 0.3];
y0 = ones(2*N, 1);
tic;
[t, y] = ode45(cartfun, tspan, y0);
toc;
c = y(:, 1:N); a = y(:, 1+N:2*N);
figure
mesh(c);
figure
mesh(a);
%% [car2rhs|car2jac] Mass matrix version of cartoon model
N = 5;
I = 0.1;
D = 0.5;
car2fun = @(t,y) car2rhs(t,y,I,D);
tspan = [0 0.3];
y0 = [ones(2*N,1); 0.5*ones(3*N,1)];
M = sparse((1:2*N), (1:2*N), ones(1,2*N), 5*N, 5*N);    % mass matrix
% M = diag([ones(1,2*N) zeros(1,3*N)]);
J = @car2jac;   % Jacobian matrix
options = odeset('Mass', M, 'Jacobian', J); % set Mass matrix and Jacobian matrix
[t, y] = ode15s(car2fun, tspan, y0, options);   % ode solver - ode15s
c = y(:, 1:N); a = y(:, 1+N:2*N);
figure
mesh(c);    % mesh figure of c
figure
mesh(a);    % mesh figure of a
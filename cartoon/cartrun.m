%CARTRUN run of cartoon model
%% [cartrhs]
% N = 5;
% I = 0.1;
% D = 0.5;
% cartfun = @(t,y) cartrhs(t,y,I,D);
% tspan = [0 0.3];
% y0 = ones(2*N, 1);
% tic;
% [t, y] = ode45(cartfun, tspan, y0);
% toc;
% c = y(:, 1:N); a = y(:, 1+N:2*N);
% figure
% mesh(c);
% figure
% mesh(a);
%% [car2rhs|car2jac] Mass matrix version of cartoon model
clear all; close all;
% symbolic expressions of Cartoon model
global c a cs V phi j % global variables of cartoon model
global y dy dF        % function element of the model
global N M I D        % parameters of the model
N = 5;      % number of spatial intervals
M = 3;      % number of radial intervals
I = 0.01;   % current
D = 0.5;    % parameter D describing the counter-ion
car2fun;    % calculation of symbolic function of cartoon model
% ode15s function used to solve DAEs and PDEs hybrid problem
odefun = @(t,y) car2rhs(t,y);
tspan = [0 0.3];
y0 = [ones((1+M)*N,1); 0.5*ones(3*N,1)];
Mas = sparse((1:(1+M)*N), (1:(1+M)*N), ones(1,(1+M)*N), (4+M)*N, (4+M)*N);    % mass matrix
% Mas = diag([ones(1,2*N) zeros(1,3*N)]);
Jac = @(t,y) car2jac(t,y);   % Jacobian matrix
options = odeset('Mass', Mas, 'Jacobian', Jac); % set Mass matrix and Jacobian matrix, 
tic
[t_, y_] = ode15s(odefun, tspan, y0, options);   % ode solver - ode15s
etime = toc
c_ = y_(:, 1:N); a_ = y_(:, 1+N:2*N);
figure
mesh(c_);    % mesh figure of c
figure
mesh(a_);    % mesh figure of a
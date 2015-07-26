% ode test
%% ode test of original func
figure;
a = 1; b = 1;
[t, y] = ode45(@(t,y) func(t,y,a,b), [1,4], 1);
subplot(2,2,1);
plot(t, y, 'linewidth', 2); grid on
title(['a = ' num2str(a) ', b = ' num2str(b)]); xlabel('t'); ylabel('y');
xlim([1 4]); ylim([0 40]);
a = 1; b = 2;
[t, y] = ode45(@(t,y) func(t,y,a,b), [1,4], 1);
subplot(2,2,(a-1)*2+b);
plot(t, y, 'linewidth', 2); grid on
title(['a = ' num2str(a) ', b = ' num2str(b)]); xlabel('t'); ylabel('y');
xlim([1 4]); ylim([0 40]);
a = 2; b = 1;
[t, y] = ode45(@(t,y) func(t,y,a,b), [1,4], 1);
subplot(2,2,(a-1)*2+b);
plot(t, y, 'linewidth', 2); grid on
title(['a = ' num2str(a) ', b = ' num2str(b)]); xlabel('t'); ylabel('y');
xlim([1 4]); ylim([0 40]);
a = 2; b = 2;
[t, y] = ode45(@(t,y) func(t,y,a,b), [1,4], 1);
subplot(2,2,(a-1)*2+b);
plot(t, y, 'linewidth', 2); grid on
title(['a = ' num2str(a) ', b = ' num2str(b)]); xlabel('t'); ylabel('y');
xlim([1 4]); ylim([0 40]);
%% ode test of multiple equations
options = odeset('RelTol', 1e-4, 'AbsTol', [1e-4 1e-4 1e-5]);
tic
[t, y] = ode45(@rigid, [0 12], [0 1 1], options);
toc
% figure
plot(t, y(:,1), '-', t, y(:,2), '-.', t, y(:,3), '.');
%% ode45 function
% global N;     % N subintervals
% global h;     % spatial intervals
% global eps;   % epsilon - dimensionless parameter
N = 20;     % N subintervals
eps = 1;    % epsilon - dimensionless parameter
tspan = [0 5];  % tspan [t0 tf]
u0 = 1*ones(N,1);   % initial value at t0

err = 1e-2; % error of tolarence

% options = odeset('RelTol', err, 'AbsTol', err*ones(1, N));
tic
% [t, y] = ode45(@pde_samp, [0 5], 1*ones(1, N), options);
[t, u] = ode45(@(t,u) pde_samp(t,u,eps), tspan, u0);
ptime = toc
% figure
% plot(t, u, 'linewidth', 2); grid on
% figure
% mesh(u); 
% xlabel('x'); ylabel('t'); zlabel('u');
%% ode15s function using Jacobian Matrix
% global N;     % N subintervals
% global h;     % spatial intervals
% global eps;   % epsilon - dimensionless parameter
N = 20;     % N subintervals
eps = 1;    % epsilon - dimensionless parameter
tspan = [0 5];  % tspan [t0 tf]
u0 = 1*ones(N,1);   % initial value at t0

err = 1e-2; % error of tolarence

% options = odeset('RelTol', err, 'AbsTol', err*ones(1, N));
options = odeset('Jacobian', @(t,u) pde_jac(t,u,eps));
tic
% [t, y] = ode45(@pde_samp, [0 5], 1*ones(1, N), options);
[t, u] = ode15s(@(t,u) pde_samp(t,u,eps), tspan, u0, options);  % Jacobian Matrix included
% [t, u] = ode15s(@(t,u) pde_samp(t,u,eps), tspan, u0);           % Jacobian Matrix not included
ptime = toc
% figure
% plot(t, u, 'linewidth', 2); grid on
% figure
% mesh(u); 
% xlabel('x'); ylabel('t'); zlabel('u');
%% Symbolic operations in MATLAB
N = 4;
u = sym('u', [N 1]);
du = sym('du', [N 1]);
du(1) = (-3*u(1)+u(2))/h^2 - (u(1)^3-u(1))/eps^2;
for i = 2:N-1
    du(i) = (u(i+1)-2*u(i)+u(i-1))/h^2 - (u(i)^3-u(i))/eps^2;
end
du(N) = (-u(N)+u(N-1))/h^2 - (u(N)^3-u(N))/eps^2
latex(du)
dR = jacobian(du, u)
latex(dR)
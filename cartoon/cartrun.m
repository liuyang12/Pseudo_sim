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
% N = 5;
% I = 0.1;
% D = 0.5;
% 
% h = 1/N;
% % symbolic expressions of Cartoon model
% global c a V phi j
% global y dy dF
% c = sym('c_', [N 1]);       % c(x,t) - Lithium electrolyte concentration
% a = sym('a_', [N 1]);       % a(x,t) - average concentration in the grain
% V = sym('V_', [N 1]);       % V(x,t) - solid phase voltage
% phi = sym('phi_', [N 1]);   % phi(x,t) - electrolyte potential
% j = sym('j_', [N 1]);       % j(x,t) - current equivalent concentration of grian lithium to electrolyte lithium (could be negative)
% 
% y = sym('y', [5*N 1]);
% dy = sym('dy', [5*N 1]);
% dF = sym('dF', [5*N 5*N]);
% % boundary conditions for time-independent variables
% Vg0 = -V(1);                         % V(0) - ghost point 0
% Vg1 = V(N);                          % V(N+1) - ghost point 1
% phig0 = phi(1);                      % phi(0) - ghost point 0
% phig1 = phi(N)+I*h/(2*c(N)+I*h/2);   % phi(N+1) - ghost point 1
% cg0 = c(1);                          % c(1) - ghost point 0
% cg1 = c(N)+I*h/2;                    % c(N+1) - ghost point 1
% % differential equations [time-dependent variables] 2*N
% % equation (5)/(6)
% dc(1) = D*(c(2)-2*c(1)+cg0)/h^2 - D*((c(2)+c(1))/2*(phi(2)-phi(1))/h - (c(1)+cg0)/2*(phi(1)-phig0)/h);
% for i = 2:N-1
%     dc(i) = D*(c(i+1)-2*c(i)+c(i-1))/h^2 - D*((c(i+1)+c(i))/2*(phi(i+1)-phi(i))/h - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))/h);
% end
% dc(N) = D*(cg1-2*c(N)+c(N-1))/h^2 - D*((cg1+c(N))/2*(phig1-phi(N))/h - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))/h);
% % equation (7)
% for i = 1:N
%     da(i) = -j(i);
% end
% dy(1:N,1) = dc';
% dy(1+N:N+N,1) = da';
% % main equations [time-independent variables] 3*N
% % equation (8)
% dy(1+2*N) = (V(2)-2*V(1)+Vg0)/h^2 - j(1);
% for i = 2:N-1
%     dy(i+2*N) = (V(i+1)-2*V(i)+V(i-1))/h^2 - j(i);
% end
% dy(N+2*N) = (Vg1-2*V(N)+V(N-1))/h^2 - j(N);
% % equation (4*)
% for i = 1:N
%     dy(i+3*N)=j(i)-sinh(V(i)-phi(i)-1+log(a(i))-log(c(i)));
% end
% % equation (6*)
% % dy(1+4*N) = (1-D)*(c(2)-2*c(1)+cg0)/h^2 + (1+D)/h*((c(2)+c(1))/2*(phi(2)-phi(1)) - (c(1)+cg0)/2*(phi(1)-phig0)) - (V(2)-2*V(1)+Vg0)/h^2;
% % for i = 2:N-1
% %     dy(i+4*N) = (1-D)*(c(i+1)-2*c(i)+c(i-1))/h^2 + (1+D)/h*((c(i+1)+c(i))/2*(phi(i+1)-phi(i)) - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))) - (V(i+1)-2*V(i)+V(i-1))/h^2;
% % end
% % dy(N+4*N) = (1-D)*(cg1-2*c(N)+c(N-1))/h^2 + (1+D)/h*((cg1+c(N))/2*(phig1-phi(N)) - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))) - (Vg1-2*V(N)+V(N-1))/h^2;
% % according to equation (8), (V(i+1)-2*V(i)+V(i-1))/h^2 can be replaced as j(i)
% dy(1+4*N) = (1-D)*(c(2)-2*c(1)+cg0)/h^2 + (1+D)/h*((c(2)+c(1))/2*(phi(2)-phi(1)) - (c(1)+cg0)/2*(phi(1)-phig0)) - j(1);
% for i = 2:N-1
%     dy(i+4*N) = (1-D)*(c(i+1)-2*c(i)+c(i-1))/h^2 + (1+D)/h*((c(i+1)+c(i))/2*(phi(i+1)-phi(i)) - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))) - j(i);
% end
% dy(N+4*N) = (1-D)*(cg1-2*c(N)+c(N-1))/h^2 + (1+D)/h*((cg1+c(N))/2*(phig1-phi(N)) - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))) - j(N);
% 
% y = [c;a;V;phi;j];
% dF = jacobian(dy, y);       % Jacobian matrix of cartoon model symbolically
clear all; close all;
% symbolic expressions of Cartoon model
global c a V phi j  % global variables of cartoon model
global y dy dF      % function element of the model
global N I D        % parameters of the model
N = 5;      % number of intervals
I = 0.1;    % current
D = 0.5;    % parameter D describing the counter-ion
car2fun;    % calculation of symbolic function of cartoon model
% ode15s function used to solve DAEs and PDEs hybrid problem
odefun = @(t,y) car2rhs(t,y);
tspan = [0 1];
y0 = [ones(2*N,1); 0.5*ones(3*N,1)];
M = sparse((1:2*N), (1:2*N), ones(1,2*N), 5*N, 5*N);    % mass matrix
% M = diag([ones(1,2*N) zeros(1,3*N)]);
J = @(t,y) car2jac(t,y);   % Jacobian matrix
options = odeset('Mass', M, 'Jacobian', J); % set Mass matrix and Jacobian matrix, 
tic
[t_, y_] = ode15s(odefun, tspan, y0, options);   % ode solver - ode15s
etime = toc
c_ = y_(:, 1:N); a_ = y_(:, 1+N:2*N);
figure
mesh(c_);    % mesh figure of c
figure
mesh(a_);    % mesh figure of a
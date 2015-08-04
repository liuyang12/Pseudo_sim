function [ dy ] = cartrhs( t, y )
%CARTRHS RHS of cartoon model by Brian Wetton, July 27, 2015
%   CARTRHS - implementation by Yang Liu, July 31, 2015
%   Detailed explanation goes here
%   y is a column vector dimensioned 2N*1 [c(1:N) a(1:N)]
%   c(x,t) - Lithium electrolyte concentration
%   a(x,t) - average concentration in the grain
%   V(x,t) - solid phase voltage
%   phi(x,t) - electrolyte potential
%   j(x,t) - current equivalent concentration of grian lithium to electrolyte lithium (could be negative)
%   q(x,t) - flux of lithium in electrolyte
%   I(x,t) - solid phase current
N = max(size(y))/2; % N subintervals for c and a.
h = 1/N;            % spatial step length 
%% initialization of time-dependent variables
c = y(1:N);         % c(x,t) - Lithium electrolyte concentration, pay attention to the order of y
a = y(N+1:2*N);     % a(x,t) - average concentration in the grain
%% caculation of time-independent variables
V0 = 0.5*ones(1,N);     % absolute initial value of vector V
phi0 = 0.5*ones(1,N);   % absolute initial value of vector phi
j0 = 0.5*ones(1,N);     % absolute initial value of vector j
x_0 = [V0 phi0 j0];     % absolute initial values of vector x
I = 2.6;    % applied current - 1C current of 2600mAh battery
D = 0.5;    % parameter D
persistent x0;          % persistant vector x0 [use the result of Newton method as initial vector of next time step]
if isempty(x0)
    x0 = x_0;           % absolute initial value of vecor x0
end
V = sym('V_', [1 N]);       % V(x,t) - solid phase voltage
phi = sym('phi_', [1 N]);   % phi(x,t) - electrolyte potential
j = sym('j_', [1 N]);       % j(x,t) - current equivalent concentration of grian lithium to electrolyte lithium (could be negative)
% syms V0 V1 phi0 phi1 c0 c1  % ghost point of V, phi and c
% boundary conditions for time-independent variables
Vg0 = -V(1);                         % V(0)
Vg1 = V(N);                          % V(N+1)
phig0 = phi(1);                      % phi(0)
phig1 = phi(N)+I*h/(2*c(N)+I*h/2);   % phi(N+1)
cg0 = c(1);                          % c(1)
cg1 = c(N)+I*h/2;                    % c(N+1)
% pay attention to the order of the equations - fx
% main equations [time-independent variables] 3*N
% equation (8)
fx(1) = (V(2)-2*V(1)+Vg0)/h^2 - j(1);
for i = 2:N-1
    fx(i) = (V(i+1)-2*V(i)+V(i-1))/h^2 - j(i);
end
fx(N) = (Vg1-2*V(N)+V(N-1))/h^2 - j(N);
% equation (4*)
for i = 1:N
    fx(i+N)=j(i)-sinh(V(i)-phi(i)-1+log(a(i))-log(c(i)));
end
% equation (6*)
% fx(1+2*N) = (1-D)*(c(2)-2*c(1)+cg0)/h^2 + (1+D)/h*((c(2)+c(1))/2*(phi(2)-phi(1)) - (c(1)+cg0)/2*(phi(1)-phig0)) - (V(2)-2*V(1)+Vg0)/h^2;
% for i = 2:N-1
%     fx(i+2*N) = (1-D)*(c(i+1)-2*c(i)+c(i-1))/h^2 + (1+D)/h*((c(i+1)+c(i))/2*(phi(i+1)-phi(i)) - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))) - (V(i+1)-2*V(i)+V(i-1))/h^2;
% end
% fx(N+2*N) = (1-D)*(cg1-2*c(N)+c(N-1))/h^2 + (1+D)/h*((cg1+c(N))/2*(phig1-phi(N)) - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))) - (Vg1-2*V(N)+V(N-1))/h^2;
% according to equation (8), (V(i+1)-2*V(i)+V(i-1))/h^2 can be replaced as j(i)
fx(1+2*N) = (1-D)*(c(2)-2*c(1)+cg0)/h^2 + (1+D)/h*((c(2)+c(1))/2*(phi(2)-phi(1)) - (c(1)+cg0)/2*(phi(1)-phig0)) - j(1);
for i = 2:N-1
    fx(i+2*N) = (1-D)*(c(i+1)-2*c(i)+c(i-1))/h^2 + (1+D)/h*((c(i+1)+c(i))/2*(phi(i+1)-phi(i)) - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))) - j(i);
end
fx(N+2*N) = (1-D)*(cg1-2*c(N)+c(N-1))/h^2 + (1+D)/h*((cg1+c(N))/2*(phig1-phi(N)) - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))) - j(N);
x = [V phi j];
[xi i] = newton(fx, x, x0); % Newton Method of solving nonlinear equations
V1 = xi(1:N);
phi1 = xi(1+N:N+N);
j1 = xi(1+2*N:N+2*N);
V1g0 = -V1(1);                         % V(0)
V1g1 = V1(N);                          % V(N+1)
phi1g0 = phi1(1);                      % phi(0)
phi1g1 = phi1(N)+I*h/(2*c(N)+I*h/2);   % phi(N+1)
c1g0 = c(1);                           % c(1)
c1g1 = c(N)+I*h/2;                     % c(N+1)
x0 = xi;    % use the result of Newton method as initial vector of next time step

% differential equations [time-dependent variables] 2*N
% equation (5)/(6)
dc(1) = D*(c(2)-2*c(1)+c1g0)/h^2 - D*((c(2)+c(1))/2*(phi1(2)-phi1(1))/h - (c(1)+c1g0)/2*(phi1(1)-phi1g0)/h);
for i = 2:N-1
    dc(i) = D*(c(i+1)-2*c(i)+c(i-1))/h^2 - D*((c(i+1)+c(i))/2*(phi1(i+1)-phi1(i))/h - (c(i)+c(i-1))/2*(phi1(i)-phi1(i-1))/h);
end
dc(N) = D*(c1g1-2*c(N)+c(N-1))/h^2 - D*((c1g1+c(N))/2*(phi1g1-phi1(N))/h - (c(N)+c(N-1))/2*(phi1(N)-phi1(N-1))/h);
% equation (7)
for i = 1:N
    da(i) = -j1(i);
end
y = [c;a];
dy = [dc';da'];
end
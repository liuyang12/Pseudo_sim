% symbolic function of Cartoon model

% symbolic expressions of Cartoon model
global c a V phi j
global y dy dF
global N I D % parameters of the model
h = 1/N;    % spatial step
c = sym('c_', [N 1]);       % c(x,t) - Lithium electrolyte concentration
a = sym('a_', [N 1]);       % a(x,t) - average concentration in the grain
V = sym('V_', [N 1]);       % V(x,t) - solid phase voltage
phi = sym('phi_', [N 1]);   % phi(x,t) - electrolyte potential
j = sym('j_', [N 1]);       % j(x,t) - current equivalent concentration of grian lithium to electrolyte lithium (could be negative)

y = sym('y', [5*N 1]);
dy = sym('dy', [5*N 1]);
dF = sym('dF', [5*N 5*N]);
% boundary conditions for time-independent variables
Vg0 = -V(1);                         % V(0) - ghost point 0
Vg1 = V(N);                          % V(N+1) - ghost point 1
phig0 = phi(1);                      % phi(0) - ghost point 0
phig1 = phi(N)+I*h/(2*c(N)+I*h/2);   % phi(N+1) - ghost point 1
cg0 = c(1);                          % c(1) - ghost point 0
cg1 = c(N)+I*h/2;                    % c(N+1) - ghost point 1
% differential equations [time-dependent variables] 2*N
% equation (5)/(6)
dc(1) = D*(c(2)-2*c(1)+cg0)/h^2 - D*((c(2)+c(1))/2*(phi(2)-phi(1))/h - (c(1)+cg0)/2*(phi(1)-phig0)/h);
for i = 2:N-1
    dc(i) = D*(c(i+1)-2*c(i)+c(i-1))/h^2 - D*((c(i+1)+c(i))/2*(phi(i+1)-phi(i))/h - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))/h);
end
dc(N) = D*(cg1-2*c(N)+c(N-1))/h^2 - D*((cg1+c(N))/2*(phig1-phi(N))/h - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))/h);
% equation (7)
for i = 1:N
    da(i) = -j(i);
end
dy(1:N,1) = dc';
dy(1+N:N+N,1) = da';
% main equations [time-independent variables] 3*N
% equation (8)
dy(1+2*N) = (V(2)-2*V(1)+Vg0)/h^2 - j(1);
for i = 2:N-1
    dy(i+2*N) = (V(i+1)-2*V(i)+V(i-1))/h^2 - j(i);
end
dy(N+2*N) = (Vg1-2*V(N)+V(N-1))/h^2 - j(N);
% equation (4*)
for i = 1:N
    dy(i+3*N)=j(i)-sinh(V(i)-phi(i)-1+log(a(i))-log(c(i)));
end
% equation (6*)
% dy(1+4*N) = (1-D)*(c(2)-2*c(1)+cg0)/h^2 + (1+D)/h*((c(2)+c(1))/2*(phi(2)-phi(1)) - (c(1)+cg0)/2*(phi(1)-phig0)) - (V(2)-2*V(1)+Vg0)/h^2;
% for i = 2:N-1
%     dy(i+4*N) = (1-D)*(c(i+1)-2*c(i)+c(i-1))/h^2 + (1+D)/h*((c(i+1)+c(i))/2*(phi(i+1)-phi(i)) - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))) - (V(i+1)-2*V(i)+V(i-1))/h^2;
% end
% dy(N+4*N) = (1-D)*(cg1-2*c(N)+c(N-1))/h^2 + (1+D)/h*((cg1+c(N))/2*(phig1-phi(N)) - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))) - (Vg1-2*V(N)+V(N-1))/h^2;
% according to equation (8), (V(i+1)-2*V(i)+V(i-1))/h^2 can be replaced as j(i)
dy(1+4*N) = (1-D)*(c(2)-2*c(1)+cg0)/h^2 + (1+D)/h*((c(2)+c(1))/2*(phi(2)-phi(1)) - (c(1)+cg0)/2*(phi(1)-phig0)) - j(1);
for i = 2:N-1
    dy(i+4*N) = (1-D)*(c(i+1)-2*c(i)+c(i-1))/h^2 + (1+D)/h*((c(i+1)+c(i))/2*(phi(i+1)-phi(i)) - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))) - j(i);
end
dy(N+4*N) = (1-D)*(cg1-2*c(N)+c(N-1))/h^2 + (1+D)/h*((cg1+c(N))/2*(phig1-phi(N)) - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))) - j(N);

y = [c;a;V;phi;j];
dF = jacobian(dy, y);       % Jacobian matrix of cartoon model symbolically
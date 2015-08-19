% symbolic function of Cartoon model

% symbolic expressions of Cartoon model
global c a cs V phi j
global y dy dF
global N M I D % parameters of the model
h = 1/N;    % spatial step
hm = 1/M;   % radial step
c = sym('c_', [N 1]);       % c(x,t) - Lithium electrolyte concentration [c_e(x,t)]
a = sym('a_', [N 1]);       % a(x,t) - average concentration in the grain [c_ss(x,t)]
cs = sym('c_s', [N M]);     % cs(x,r,t) - solid phase concentration with radial distribution [c_s(x,t)]
                            %             where c_ss(x,t) = c_s(x,Rp,t)
V = sym('V_', [N 1]);       % V(x,t) - solid phase voltage
phi = sym('phi_', [N 1]);   % phi(x,t) - electrolyte potential
j = sym('j_', [N 1]);       % j(x,t) - current equivalent concentration of grian lithium to electrolyte lithium (could be negative)

y = sym('y', [(4+M)*N 1]);
dy = sym('dy', [(4+M)*N 1]);
dF = sym('dF', [(4+M)*N (4+M)*N]);
% boundary conditions for time-independent variables
Vg0 = -V(1);                         % V(0) - ghost point 0
Vg1 = V(N);                          % V(N+1) - ghost point 1
phig0 = phi(1);                      % phi(0) - ghost point 0
phig1 = phi(N)+I*h/(2*c(N)+I*h/2);   % phi(N+1) - ghost point 1
cg0 = c(1);                          % c(1) - ghost point 0
cg1 = c(N)+I*h/2;                    % c(N+1) - ghost point 1
csg0 = cs(:,1);                     % cs(0) - ghost vector 0
csg1 = cs(:,M) - hm*j;              % cs(N+1) - ghost vector 1
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
% equation (1.4) in [Ref. 3]
dcs(:,1) = 1/(4*1^2*hm)*((2*1+1)^2*(cs(:,1+1)-cs(:,1)) - (2*1-1)^2*(cs(:,1)-csg0));
for m = 2:M-1
   dcs(:,m) = 1/(4*m^2*hm)*((2*m+1)^2*(cs(:,m+1)-cs(:,m)) - (2*m-1)^2*(cs(:,m)-cs(:,m-1)));
end
dcs(:,M) = 1/(4*M^2*hm)*((2*M+1)^2*(csg1-cs(:,M)) - (2*M-1)^2*(cs(:,M)-cs(:,M-1)));

dy(1:N,1) = dc';
dy(1+N:M*N+N,1) = dcs(:)';
% main equations [time-independent variables] 3*N
% equation (8)
dy(1+(1+M)*N) = (V(2)-2*V(1)+Vg0)/h^2 - j(1);
for i = 2:N-1
    dy(i+(1+M)*N) = (V(i+1)-2*V(i)+V(i-1))/h^2 - j(i);
end
dy(N+(1+M)*N) = (Vg1-2*V(N)+V(N-1))/h^2 - j(N);
% equation (4*)
a = 1/2*(cs(:,M) + csg1);   % css(x,t) = cs(x,Rp,t)
for i = 1:N
    dy(i+(2+M)*N)=j(i)-sinh(V(i)-phi(i)-1+log(a(i))-log(c(i)));
end
% equation (6*)
% dy(1+4*N) = (1-D)*(c(2)-2*c(1)+cg0)/h^2 + (1+D)/h*((c(2)+c(1))/2*(phi(2)-phi(1)) - (c(1)+cg0)/2*(phi(1)-phig0)) - (V(2)-2*V(1)+Vg0)/h^2;
% for i = 2:N-1
%     dy(i+4*N) = (1-D)*(c(i+1)-2*c(i)+c(i-1))/h^2 + (1+D)/h*((c(i+1)+c(i))/2*(phi(i+1)-phi(i)) - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))) - (V(i+1)-2*V(i)+V(i-1))/h^2;
% end
% dy(N+4*N) = (1-D)*(cg1-2*c(N)+c(N-1))/h^2 + (1+D)/h*((cg1+c(N))/2*(phig1-phi(N)) - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))) - (Vg1-2*V(N)+V(N-1))/h^2;
% according to equation (8), (V(i+1)-2*V(i)+V(i-1))/h^2 can be replaced as j(i)
dy(1+(3+M)*N) = (1-D)*(c(2)-2*c(1)+cg0)/h^2 + (1+D)/h*((c(2)+c(1))/2*(phi(2)-phi(1)) - (c(1)+cg0)/2*(phi(1)-phig0)) - j(1);
for i = 2:N-1
    dy(i+(3+M)*N) = (1-D)*(c(i+1)-2*c(i)+c(i-1))/h^2 + (1+D)/h*((c(i+1)+c(i))/2*(phi(i+1)-phi(i)) - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))) - j(i);
end
dy(N+(3+M)*N) = (1-D)*(cg1-2*c(N)+c(N-1))/h^2 + (1+D)/h*((cg1+c(N))/2*(phig1-phi(N)) - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))) - j(N);

y = [c;cs(:);V;phi;j];
dF = jacobian(dy, y);       % Jacobian matrix of cartoon model symbolically
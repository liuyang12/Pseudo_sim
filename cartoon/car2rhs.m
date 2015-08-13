function [ dydt ] = car2rhs( t, y1 )
%CAR2RHS Mass matrix version of cartoon model
%   CAR2RHS - vector y contains all the 5*N (c,a,V,phi,j) variables in the model
%   
%   I = 2.6;    % applied current - 1C current of 2600mAh battery
%   D = 0.5;    % parameter D
global c a V phi j
global y dy dF

% N = max(size(y1))/5;
% h = 1/N;
% % pay attention to the order of the variables, especially the time-dependent and
% % time-dependent variables
% % c = y(1:N);
% % a = y(1+N:N+N);
% % V = y(1+2*N:N+2*N);
% % phi = y(1+3*N:N+3*N);
% % j = y(1+4*N:N+4*N);
% 
% c = sym('c_', [N 1]);       % c(x,t) - Lithium electrolyte concentration
% a = sym('a_', [N 1]);       % a(x,t) - average concentration in the grain
% V = sym('V_', [N 1]);       % V(x,t) - solid phase voltage
% phi = sym('phi_', [N 1]);   % phi(x,t) - electrolyte potential
% j = sym('j_', [N 1]);       % j(x,t) - current equivalent concentration of grian lithium to electrolyte lithium (could be negative)
% 
% % persistent dy;  % use persistent/global symbolic function of dy
% % if isempty(dy)
%     % boundary conditions for time-independent variables
%     Vg0 = -V(1);                         % V(0) - ghost point 0
%     Vg1 = V(N);                          % V(N+1) - ghost point 1
%     phig0 = phi(1);                      % phi(0) - ghost point 0
%     phig1 = phi(N)+I*h/(2*c(N)+I*h/2);   % phi(N+1) - ghost point 1
%     cg0 = c(1);                          % c(1) - ghost point 0
%     cg1 = c(N)+I*h/2;                    % c(N+1) - ghost point 1
%     % differential equations [time-dependent variables] 2*N
%     % equation (5)/(6)
%     dc(1) = D*(c(2)-2*c(1)+cg0)/h^2 - D*((c(2)+c(1))/2*(phi(2)-phi(1))/h - (c(1)+cg0)/2*(phi(1)-phig0)/h);
%     for i = 2:N-1
%         dc(i) = D*(c(i+1)-2*c(i)+c(i-1))/h^2 - D*((c(i+1)+c(i))/2*(phi(i+1)-phi(i))/h - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))/h);
%     end
%     dc(N) = D*(cg1-2*c(N)+c(N-1))/h^2 - D*((cg1+c(N))/2*(phig1-phi(N))/h - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))/h);
%     % equation (7)
%     for i = 1:N
%         da(i) = -j(i);
%     end
%     dy(1:N,1) = dc';
%     dy(1+N:N+N,1) = da';
%     % main equations [time-independent variables] 3*N
%     % equation (8)
%     dy(1+2*N) = (V(2)-2*V(1)+Vg0)/h^2 - j(1);
%     for i = 2:N-1
%         dy(i+2*N) = (V(i+1)-2*V(i)+V(i-1))/h^2 - j(i);
%     end
%     dy(N+2*N) = (Vg1-2*V(N)+V(N-1))/h^2 - j(N);
%     % equation (4*)
%     for i = 1:N
%         dy(i+3*N)=j(i)-sinh(V(i)-phi(i)-1+log(a(i))-log(c(i)));
%     end
%     % equation (6*)
%     % dy(1+4*N) = (1-D)*(c(2)-2*c(1)+cg0)/h^2 + (1+D)/h*((c(2)+c(1))/2*(phi(2)-phi(1)) - (c(1)+cg0)/2*(phi(1)-phig0)) - (V(2)-2*V(1)+Vg0)/h^2;
%     % for i = 2:N-1
%     %     dy(i+4*N) = (1-D)*(c(i+1)-2*c(i)+c(i-1))/h^2 + (1+D)/h*((c(i+1)+c(i))/2*(phi(i+1)-phi(i)) - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))) - (V(i+1)-2*V(i)+V(i-1))/h^2;
%     % end
%     % dy(N+4*N) = (1-D)*(cg1-2*c(N)+c(N-1))/h^2 + (1+D)/h*((cg1+c(N))/2*(phig1-phi(N)) - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))) - (Vg1-2*V(N)+V(N-1))/h^2;
%     % according to equation (8), (V(i+1)-2*V(i)+V(i-1))/h^2 can be replaced as j(i)
%     dy(1+4*N) = (1-D)*(c(2)-2*c(1)+cg0)/h^2 + (1+D)/h*((c(2)+c(1))/2*(phi(2)-phi(1)) - (c(1)+cg0)/2*(phi(1)-phig0)) - j(1);
%     for i = 2:N-1
%         dy(i+4*N) = (1-D)*(c(i+1)-2*c(i)+c(i-1))/h^2 + (1+D)/h*((c(i+1)+c(i))/2*(phi(i+1)-phi(i)) - (c(i)+c(i-1))/2*(phi(i)-phi(i-1))) - j(i);
%     end
%     dy(N+4*N) = (1-D)*(cg1-2*c(N)+c(N-1))/h^2 + (1+D)/h*((cg1+c(N))/2*(phig1-phi(N)) - (c(N)+c(N-1))/2*(phi(N)-phi(N-1))) - j(N);
% % end
% 
% % symbolical to numberical
% y = [c;a;V;phi;j];
dydt = double(subs(dy, y, y1)); % RHS of cartoon model numerically
end


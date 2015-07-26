function du = pde_samp(t, u, eps)
% pde function of example problem
N = max(size(u));   % N subintervals
h = 1/N;            % spatial intervals
% eps = 1;            % epsilon - dimensionless parameter
% du = zeros(N, 1);
% du(1) = (-3*u(1)+u(2))/h^2 - (u(1)^3-u(1))/eps^2;
% for i = 2:N-1
%     du(i) = (u(i+1)-2*u(i)+u(i-1))/h^2 - (u(i)^3-u(i))/eps^2;
% end
% du(N) = (-u(N)+u(N-1))/h^2 - (u(N)^3-u(N))/eps^2;
%% another implementation of RHS of Allen-Cahn equation by Brian Wetton
% useful structure in diiferencial problems
ur = zeros(N,1);
ur(1:N-1) = u(2:N);
ur(N) = u(N);   % boundary condition 2

ul = zeros(N,1);
ul(2:N) = u(1:N-1);
ul(1) = -u(1);  % boundary condition 1

du = (ul+ur-2*u)/h^2 - (u.^3 - u)/eps^2;  % RHS of PDE

end

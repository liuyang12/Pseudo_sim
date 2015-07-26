function [ J ] = pde_jac(t, u, eps)
% jacobian matrix of demo pde function
N = max(size(u));   % N subintervals
h = 1/N;            % spatial intervals
% eps = 1;            % epsilon - dimensionless parameter

i(1:N-1) = 1:N-1;       
j(1:N-1) = 2:N;
v(1:N-1) = 1/h^2;
i(N:2*N-2) = 2:N;       
j(N:2*N-2) = 1:N-1;
v(N:2*N-2) = 1/h^2;
i(2*N-1:3*N-2) = 1:N;   
j(2*N-1:3*N-2) = 1:N;
v(2*N-1) = -3/h^2;
v(2*N:3*N-3) = -2/h^2;
v(3*N-2) = -1/h^2;
v(2*N-1:3*N-2) = v(2*N-1:3*N-2) - (3*(u.^2)' - 1)/eps^2;
J = sparse(i, j, v, N, N);

% J = sparse(N, N);   % sparase matrix 
% for j = 1:N     % calculation of Jacobian Matrix for each row_j
%     if j == 1
%         J(j, j) = -3/h^2;
%         J(j, j+1) = 1/h^2;
%     elseif j == N
%         J(j, j) = -1/h^2;
%         J(j, j-1) = 1/h^2;
%     else
%         J(j, j) = -2/h^2;
%         J(j, j-1) = 1/h^2;
%         J(j, j+1) = 1/h^2;
%     end
%     J(j, j) = J(j, j) - (3*u(j)^2 - 1)/eps^2;
% end

end
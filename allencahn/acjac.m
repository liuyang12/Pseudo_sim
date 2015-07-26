function [ J ] = acjac( t,u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = max(size(u));
h = 1/N;

J = sparse(N,N);

for j=1:N
    if j==1 
        J(1,1) = -3/h^2;
        J(1,2) = 1/h^2;
    elseif j==N
        J(N,N) = -1/h^2;
        J(N,N-1) = 1/h^2;
    else 
        J(j,j) = -2/h^2;
        J(j,j+1) = 1/h^2;
        J(j,j-1) = 1/h^2;
    end
    
    J(j,j) = J(j,j) -3*u(j)^2 +1;
end

end


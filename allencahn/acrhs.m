function [ fu ] = acrhs( t, u )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



N = max(size(u));
h = 1/N;

ur = zeros(N,1);
ur(1:N-1) = u(2:N);
ur(N) = u(N);

ul = zeros(N,1);
ul(2:N) = u(1:N-1);
ul(1) = -u(1);

fu = (ul+ur-2*u)/h^2 - (u.^3 - u);

end

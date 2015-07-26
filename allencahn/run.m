N = 50;

S = sparse(N,N);
for j=1:N
    S(j,j) = 1;
    
    if j ~= 1
        S(j,j-1) = 1;
    end
    
    if j ~= N
        S(j, j+1) = 1;
    end
end

options = odeset('Jacobian',@acjac,'Jpattern',S);

tic
[t y] = ode15s(@acrhs, [0 5], ones(N,1), options);
toc
temp = size(y);

disp(sprintf('Time steps: %d \n',temp(1)));
plot(t(2:end), t(2:end)-t(1:end-1))
title('Time Steps')
xlabel('t')

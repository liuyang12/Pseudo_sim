function dy = rigid(t, y)
% demo function of multiple equations using ode45
dy = zeros(3, 1);
dy(1) = y(2) * y(3);
dy(2) = -y(1) * y(3);
dy(3) = -.51*y(1)*y(2);
end
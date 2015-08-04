function [xi,i]=newton(fx,x,x0,eps,N)
% NEWTON Newton's method to solve nonlinear equations.
%
% NEWTON(fx,x,x0,eps,N) using newton iteration to search solution 
% of nonlinear equations 'fx' in function of 'x' nearby inital values 'x0',
% with tolerance 'eps' or maximum step 'N'.
% NEWTON(fx,x,x0,eps) is the same excpet defautly set N=100.
% NEWTON(fx,x,x0) set eps=1.0e-3 and N=100.
%
% ARGUMENTS:
%       fx       The symsbolic row array of eqs.
%       x        The symsbolic variables of eqs.
%       x0       Guessed initial values of x.
%       eps      The tolorence for iteration.
%       N        The maximum steps for iteration.
% 
% Example 1:
% % Fixed points for Quadratic Integerate-and-Fire model (1-dimension)
% syms v
% fx=v*(v-5);
% x=v;
% [X1,t]=newton(fx,x,2.4,1.0e-10,100)
% [X2,t]=newton(fx,x,2.6,1.0e-10,100)
%
% Example 2:
% % Fixed point for Classic Hodgkin-Huxley model (4-dimensions)
% syms v m h n
% gNa=120.0; gK=36.0; gL=0.3;
% vNa=50; vK=-77; vL=-54.4;
% f1=-gNa*m^3*h*(v-vNa) -gK*n^4*(v-vK) -gL*(v-vL) ;
% f2=0.1*(v+40)/(1-exp(-(40+v)/10))*(1-m)  -4.0*exp(-(65+v)/18)*m;
% f3=0.07*exp(-(65+v)/20)*(1-h)  -1/(exp(-(35+v)/10)+1)*h;
% f4=0.01*(v+55)/(1-exp(-(v+55)/10))*(1-n)-0.125*exp(-(v+65)/18)*n;
% x=[v,m,h,n];
% fx=[f1,f2,f3,f4];
% [X,n]=newton(fx,x,[-100,1.0,1.0,1.0],1.0e-15,100);
% format long
% disp(n)
% disp(X)
% 
% ALGORITHM:
%        X(n+1) = X(n) - F(X(n)) / DF(X(n)) 
%                     = X(n) - inv(DF(X(n))) * F(X(n))
% where DF(Xn) is the differential value at X(n).
% The iteration stop when max(X(n+1)-X(n))<eps
% or maxmum step N reaching.
%
% AUTHOR: felonwan@gmail.com
% CREATED: 2013-09-04
% LAST MODIFIED: 2013-09-04

if nargin==3
	eps=1.0e-3;
	N=10;
elseif nargin==4
	N=10;
elseif nargin~=5
    error('Too few or too many arguments! (3--5)')
end
if isrow(x0)
    x0=x0.';
end
if isrow(x)
    x=x.';
end
if isrow(fx)
    fx=fx.';
end
dfx=jacobian(fx,x);
for i=1:N
    xi=x0-double(subs(dfx,x,x0))\double(subs(fx,x,x0)); % notice: subs give symbolic output in MATLAB
    if max(abs(xi-x0))<eps
        break;
    end
    x0=xi;
end

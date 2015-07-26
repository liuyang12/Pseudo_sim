# Pseudo_sim
simulation of Lithium ion batteries using Pseudo two-dimensional (P2D) model

# Solving PDE problem with ODE functions in MATLAB<sup>?</sup>
PDE is [Partial Differential Equation](https://en.wikipedia.org/wiki/Partial_differential_equation "PDE Wiki") for short. Similarly, ODE is  [Ordinary Differential Equation](https://en.wikipedia.org/wiki/ordinary_differential_equation "ODE Wiki") for short. There are several ODE functions in MATLAB<sup>?</sup>[^matlab] using numerical methods to approximate the result of ODE problems.
  [^matlab]: MATLAB<sup>?</sup> is a registered trademark of [MathWorks](http://www.mathworks.com "MathWorks official website")<sup>?</sup>.

## ODE functions in MATLAB
[ODE functions](http://www.mathworks.com/help/releases/R2015a/matlab/ref/ode45.html "ode functions in MATLAB document") in MATLAB are list as follow:

| solver | Problem type | Order of Accuracy | When to use |
| :-----: | :-----------: | :----------------: | :---------- |
| `ode45` | Nonstiff    | Medium            | Most of the time. This should be the first solver you try. |
| `ode15s` | stiff       | Low to medium     | If `ode45` is slow because the problem is stiff. |

<center>Table 1. ODE function list in MATLAB</center>

## Usage of `ode45` and `ode15s`
As is illustrated in MATLAB `help` and documentations, [`ode45`](http://www.mathworks.com/help/releases/R2015a/matlab/ref/ode45.html "ode45 functions in MATLAB document") and [ `ode15s`](http://www.mathworks.com/help/releases/R2015a/matlab/ref/ode15s.html "ode15s functions in MATLAB document") are similar in usage and basic structure, while the benefits and efficiency may differ when applied to different problems. The following is taking `ode45` for instance, followed by comparison between `ode45` and `ode15s`.
### Basic structure of `ode45` 
Basic structure of `ode45` is 
```matlab
% basic structure of ode45 in MATLAB
[TOUT,YOUT] = ode45(ODEFUN,TSPAN,Y0);
```
An example is as follow:
```matlab
[t, y] = ode45(@func, [0 4], 1);
```
where `func` refer to a differential function 
\begin{equation}
\frac{\mathrm{d}y}{\mathrm{d}t}=\frac{y}{t}+1
\label{equd1}
\end{equation}
as is shown in MATLAB script
```matlab
function dy = func(t, y)
% demo function of ode45
dy = y/t + 1; % function y with regard to t
end
```
Then, we plot this figure out.
<center>![demo function of ode45](http://img.blog.csdn.net/20150722123313560)</center>
<center>Figure 1. Demo function of `ode45`</center>

Let we think of some situation that parameters are given to `ODEFUN`. We take `func` for instance.
Consider the differential function \ref{equd1} with two parameters $a$ and $b$,
\begin{equation}
\frac{\mathrm{d}y}{\mathrm{d}t}=a\frac{y}{t}+b
\label{equd2}
\end{equation}
the function in MATLAB script is supposed to be
```matlab
function dy = func(t, y, a, b)
% demo function of ode45
dy = a*y/t + b; % function y with regard to t
end
```
which calls for adapt to this function with parameters.
### Parameter of `ODEFUN` in `ode45`
Let us look into `ODEFUN` in `ode45` function, we can see the original usage of `ODEFUN` goes to `@(t,y) func(t,y)` which is shorted by `@func` with default no extra parameters. And in such case that parameters are demanded, we can use 
```matlab
[t, y] = ode45(@(t,y) func(t,y,a,b), [0 4], 1);
```
with `a` and `b` defined ahead.
And we can do parameter sweeping to see how parameters $a$ and $b$ influence the distribution of the differential function \ref{equd2}.
<center>![parameter sweep of a and b](http://img.blog.csdn.net/20150722131616416)</center>
<center>Figure 2. Parameter sweep of $a$ and $b$</center>

### Nonstiff and stiff problem
As is shown in **Table 1**, `ode45` is applied to solve nonstiff problems, while `ode15s` is applied to solve stiff problems. So what is a stiff problem?
>In mathematics, a stiff equation is a differential equation for which certain numerical methods for solving the equation are numerically unstable, unless the step size is taken to be extremely small. It has proven difficult to formulate a precise definition of stiffness, but the main idea is that the equation includes some terms that can lead to rapid variation in the solution. -- [[Wikipedia *Stiff equation*](https://en.wikipedia.org/wiki/Stiff_equation)]

Thus we can see a nonstiff problem is the one that is not stiff. The following examples which is based on [MATLAB documentation](http://radio.feld.cvut.cz/matlab/techdoc/math_anal/ch_8_od8.html) from [Department of Radio Engineering, Czech Technical University](http://radio.feld.cvut.cz/) will show the characteristic of nonstiff and stiff problems.

### Nonstiff problem using `ode45`
A typical example in MATLAB is `rigidode`, which is first illustrated in Shampine and Gordon's book [*Computer Solution of Ordinary Differential Equations, 1975*](http://trove.nla.gov.au/version/12693184) solves the Euler equaions of a rigid body without external forces.
\begin{equation}
    \begin{cases}
    y_1'=y_2y_3, \\\\
    y_2'=-y_1y_3, \\\\
    y_3'=-0.51y_1y_2.
    \end{cases}
\label{equrigid}
\end{equation}
The solution of the ploblem \ref{equrigid} is shown in `rigidode` in MATLAB. The simplified implement is as follow:
```matlab
tspan = [0 12];
y0 = [0; 1; 1];
% solve the problem using ODE45
figure;
ode45(@f,tspan,y0);

function dydt = f(t,y)
dydt = [      y(2)*y(3)
             -y(1)*y(3)
        -0.51*y(1)*y(2) ];
```
As is shown in this MATLAB script, `ode45` is used to solve this equation numerically. If we just type `rigidode` in MATLAB command line, the following figure will appear.
<center>![rigidode in MATLAB using `ode45`](http://img.blog.csdn.net/20150724124937867)</center>
<center>Figure 3. `rigidode` in MATLAB using `ode45`</center>

### Stiff problem using `ode15s`
Van der Pol equation is a typical stiff problem. The differential equations are as follow
\begin{equation}
    \begin{cases}
        y_1'=y_2,\\\\
        y_2'=\mu(1-y_1^2)y_2-y_1.
    \end{cases}
\label{equvdp}
\end{equation}
which involves a constant parameter $\mu$.
And it is reformed from the second-order nonlinear ODE:
\begin{equation}
    y''-\mu (1-y^2)y'+y=0
\end{equation}
Then we analysis the MATLAB script using `ode15s`. 
```matlab
function  vdpode(MU)
if nargin < 1
   MU = 1000;     % default
end
tspan = [0; max(20,3*MU)];              % several periods
y0 = [2; 0];
options = odeset('Jacobian',@J);

[t,y] = ode15s(@f,tspan,y0,options);
figure;
plot(t,y(:,1));
title(['Solution of van der Pol Equation, \mu = ' num2str(MU)]);
xlabel('time t');
ylabel('solution y_1');
axis([tspan(1) tspan(end) -2.5 2.5]);

% Nested functions -- MU is provided by the outer function.
   function dydt = f(t,y)
      % Derivative function. MU is provided by the outer function.
      dydt = [            y(2)
         MU*(1-y(1)^2)*y(2)-y(1) ];
   end
   function dfdy = J(t,y)
      % Jacobian function. MU is provided by the outer function.
      dfdy = [         0                  1
         -2*MU*y(1)*y(2)-1    MU*(1-y(1)^2) ];
   end
end  % vdpode
```
Here we can see Sensitivity Matrix or Jacobian Matrix is applied to speed the process.
\begin{equation}
    \mathbf{f}(\mathbf{y})=[f_1(\mathbf{y}), f_2(\mathbf{y})]^{\mathrm T}, \; \text{where }\mathbf{y}=[y_1, y_2]^{\mathrm T}.
\label{equjac}
\end{equation}
hence, Jacobian Matrix is 
\begin{equation}
    \mathbb{J}=\frac{\partial \mathbf{f}}{\partial \mathbf{y}}=
    \left[
    {\begin{array}
    {{c}}{\frac{\partial f_1}{\partial y_1}}&{\frac{\partial f_1}{\partial y_2}}\\\\
    {\frac{\partial f_2}{\partial y_1}}&{\frac{\partial f_2}{\partial y_2}}
    \end{array}} 
    \right]
\label{equjacdef}
\end{equation}
Here the Jacobian Matrix is 
\begin{equation}
    \mathbb{J}=    
    \left[
    {\begin{array}
    {{c}}{0}&{1}\\\\
    {-2\mu y_1y_2-1}&{\mu (1-y_1^2)}
    \end{array}} 
    \right]
\end{equation}
Besides we can see usage of `ode15s` with Jacobain Matrix. Function `odeset` is used to include Jacobian Matrix in ODE functions.
```matlab
options = odeset('Jacobian', @J);
```
is an example.

Allen-Cahn equation $N=4$ 
\begin{equation}
  \mathbf{R}=\frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}=
  \left[\begin{array}{ccccc} \frac{{u_1} - u_1^3}{\varepsilon^2} - \frac{3{u_1} - {u_2}}{h^2} & \frac{{u_2} - u_2^3}{\varepsilon^2} + \frac{{u_1} - 2{u_2} + {u_3}}{h^2} & \frac{{u_3} - u_3^3}{\varepsilon^2} + \frac{{u_2} - 2{u_3} + {u_4}}{h^2} & \frac{{u_4} - u_4^3}{\varepsilon^2} + \frac{{u_3} - {u_4}}{h^2} \end{array}\right]^\mathrm{T}
\end{equation}

\begin{equation}
\mathbb{J}=\left(\frac{\partial R_i}{\partial u_j}\right)_{N\times N}=
  \left(\begin{array}{cccc}  - \frac{3u_1^2 - 1}{\varepsilon^2} - \frac{3}{h^2} & \frac{1}{h^2} & 0 & 0\\\\
   \frac{1}{h^2} &  - \frac{3u_2^2 - 1}{\varepsilon^2} - \frac{2}{h^2} & \frac{1}{h^2} & 0\\\\
    0 & \frac{1}{h^2} &  - \frac{3u_3^2 - 1}{\varepsilon^2} - \frac{2}{h^2} & \frac{1}{h^2}\\\\
     0 & 0 & \frac{1}{h^2} &  - \frac{3u_4^2 - 1}{\varepsilon^2} - \frac{1}{h^2}\\\\
      0 & 0 & 0 & \frac{1}{h^2} \end{array}\right)
\end{equation}
%% Various Parameters
Prob1.K.C = 1;
lal= 10e-6;
lp = 80e-6;
ls = 25e-6;
ln = 88e-6;
lco= 10e-6;
lt = lp+ls+ln;

dxp= lp/(Np+1);
dxs= ls/(Ns+1);
dxn= ln/(Nn+1);
dxal = lal/(Nal+1);
dxco = lco/(Nco+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_p = 0.385;    % Porosity on the positive electrode side
eps_s = 0.724;
eps_n = 0.485;
eps_i = [eps_p;eps_s;eps_n];

eps_fi= [0.025;0;0.0326];
brugg_p = 4;      % Bruggeman coefficient
brugg_s = 4;
brugg_n = 4;
brugg_i = [brugg_p;brugg_s;brugg_n];

Dps     = 1e-14;
Dns     = 3.9e-14;

a_p    = 885000;   % Particle surface area
a_s    = 0;
a_n    = 723600;

tplus = 0.364;    % Transference number - Not available for positive/negative electrode

k_p    = 2.334e-11;
k_s    = 0;
k_n    = 5.031e-11;

cs_maxp= 51554;
cs_maxs= 0;
cs_maxn= 30555;
cs_initp= 0.4955*51554;
cs_inits= 0;
cs_initn= 0.8551*30555;
c0     = 1000;

F      = 96487;
R      = 8.314;
%--------------------------------------------------------------------------
Tref   = 298.15;
%--------------------------------------------------------------------------
Rp     = 2e-6;

sig    = [100;0;100];

sig_eff = sig.*(1 - eps_i - eps_fi)/Prob1.K.C;

EaDs  = 5000;
Eaki  = 5000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ip     = 2:Np+1;
Is     = Np+3:Np+Ns+2;
In     = Np+Ns+4:Np+Ns+Nn+3;
Ico    = 1:Nco-1;

eps_i = [eps_p;eps_s;eps_n];
a_i   = [a_p;a_s;a_n];
cs_max= [cs_maxp;cs_maxs;cs_maxn];

%%%%%%%%%%%%%%%%%%%%
% Create Matrix A!
%%%%%%%%%%%%%%%%%%%%


i = [2:Np+1 2:Np+1 2:Np+1];
j = [1:Np 2:Np+1 3:Np+2];
s = [-ones(1,Np) ones(1,Np) zeros(1,Np)];
Ap1 = sparse(i,j,s,Np+2,Np+2);
Ap1(1,1:3) = [1 0 0];
Ap1(Np+2,Np:Np+2) = [1 -4 3];
Api1 = -dxp/sig_eff(1)*sparse(eye(Np,Np));
Api1 = [zeros(1,Np);Api1;zeros(1,Np)];
Api1 = [zeros(Np+2,1) Api1 zeros(Np+2,1)];
Api1(Np+2,Np+2) = -2*dxp/sig_eff(1);
Z = sparse(zeros(Np+2,Np+2));
si  = [-ones(1,Np) ones(1,Np) zeros(1,Np)];
Api = sparse(i,j,si,Np+2,Np+2);
Api(1,:) = [1 zeros(1,Np+1)];
Api(Np+2,end-2:end) = [1 -4 3];

Apphi1 = [Ap1 Api1;Z Api];
Prob1.M.iApphi1 = inv(Apphi1);

i = [2:Nn+1 2:Nn+1 2:Nn+1];
j = [1:Nn 2:Nn+1 3:Nn+2];
s = [-ones(1,Nn) ones(1,Nn) zeros(1,Nn)];
An1 = sparse(i,j,s,Nn+2,Nn+2);
An1(1,1:3) = [1 0 0];
An1(Nn+2,Nn:Nn+2) = [1 -4 3];
Ani1 = -dxn/sig_eff(3)*sparse(eye(Nn,Nn));
Ani1 = [zeros(1,Nn);Ani1;zeros(1,Nn)];
Ani1 = [zeros(Nn+2,1) Ani1 zeros(Nn+2,1)];
Ani1(Nn+2,Nn+2) = -2*dxn/sig_eff(3);
Z = sparse(zeros(Nn+2,Nn+2));
si  = [-ones(1,Nn) ones(1,Nn) zeros(1,Nn)];
Ani = sparse(i,j,si,Nn+2,Nn+2);
Ani(1,:) = [1 zeros(1,Nn+1)];
Ani(Nn+2,end-2:end) = [1 -4 3];

Anphi1 = [An1 Ani1;Z Ani];
Prob1.M.iAnphi1 = inv(Anphi1);

%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%
Prob1.M.fcp = zeros(1,Np+2);
Prob1.M.fcn = zeros(1,Nn+2);

%%%%%%%%%%%%

ip2 = [2:Nn+Ns+Np+3 2:Nn+Ns+Np+3 2:Nn+Ns+Np+3];
jp2 = [1:Nn+Ns+Np+2 2:Nn+Ns+Np+3 3:Nn+Ns+Np+4];
sp2 = [-ones(1,Nn+Ns+Np+2) ones(1,Nn+Ns+Np+2) zeros(1,Nn+Ns+Np+2)];
Ap2 = sparse(ip2,jp2,sp2,Np+Ns+Nn+4,Np+Ns+Nn+4);
Ap2(1,1:3) = [1 0 0];
Ap2(Np+Ns+Nn+4,end-2:end) = [1 -4 3];

Prob1.M.gp(1,1) = 0;
Prob1.M.gp(Nn+2,1) = 0;
Prob1.M.gp(Nn+Ns+3,1)= 0;
Prob1.M.gp(Np+Ns+Nn+4,1) = 0;


Ini = Np+Ns+Nn+3:-1:Np+Ns+4;
Isi = Np+Ns+2:-1:Np+3;
Ipi = Np+1:-1:2;
%%%%%%%%%%%
%Temperature parameters
%%%%%%%%%%%%%

% Temperature parameters
lambda_al = 237;
lambda_p  = 2.1;
lambda_s  = 0.16;
lambda_n  = 1.7;
lambda_co = 401;

rho_al = 2700;
rho_p  = 2500;
rho_s  = 1100;
rho_n  = 2500;
rho_co = 8940;

Cpal   = 897;
Cpp    = 700;
Cps    = 700;
Cpn    = 700;
Cpco   = 385;

sig_al = 3.55e7;
sig_co = 5.96e7;

%%%%%%%%%%%%%
% positive collector
%%%%%%%%%%%%%

iTal = [2:Nal 2:Nal 2:Nal];
jTal = [1:Nal-1 2:Nal 3:Nal+1];
sTal = [-lambda_al/2/dxal^2*ones(1,Nal-1) (rho_al*Cpal/dt+lambda_al/dxal^2)*ones(1,Nal-1) -lambda_al/2/dxal^2*ones(1,Nal-1)];
ATal = sparse(iTal,jTal,sTal,Nal+1,Nal+1);
ATal(1,1:3)  = [-3 4 -1];
ATal_full = [ATal sparse(Nal+1,Np+Ns+Nn+Nco+5)];
ATal_full(Nal+1,Nal:Nal+2) = [-lambda_al/2/dxal^2 (rho_al*Cpal/dt+lambda_al/dxal^2) -lambda_al/2/dxal^2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iTp = [2:Np+1 2:Np+1 2:Np+1];
jTp = [1:Np 2:Np+1 3:Np+2];
sTp = [-lambda_p/2/dxp^2*ones(1,Np) (rho_p*Cpp/dt+lambda_p/dxp^2)*ones(1,Np) -lambda_p/2/dxp^2*ones(1,Np)];
ATp = sparse(iTp,jTp,sTp,Np+2,Np+2);
ATp_full = [sparse(Np+2,Nal+1) ATp sparse(Np+2,Ns+Nn+Nco+3)];
ATp_full(1,Nal:Nal+4)  = [lambda_al/dxal -4*lambda_al/dxal 3*(lambda_al/dxal+lambda_p/dxp) -4*lambda_p/dxp lambda_p/dxp];
ATp_full(Np+2,Nal+Np+1:Nal+Np+5) = [lambda_p/dxp -4*lambda_p/dxp 3*(lambda_p/dxp+lambda_s/dxs) -4*lambda_s/dxs lambda_s/dxs];

iTs = [2:Ns-1 2:Ns-1 2:Ns-1];
jTs = [1:Ns-2 2:Ns-1 3:Ns];
sTs = [-lambda_s/2/dxs^2*ones(1,Ns-2) (rho_s*Cps/dt+lambda_s/dxs^2)*ones(1,Ns-2) -lambda_s/2/dxs^2*ones(1,Ns-2)];
ATs = sparse(iTs,jTs,sTs,Ns,Ns);
ATs_full = [sparse(Ns,Nal+Np+3) ATs sparse(Ns,Nn+Nco+3)];
ATs_full(1,Nal+Np+3:Nal+Np+5) = [-lambda_s/2/dxs^2 (rho_s*Cps/dt+lambda_s/dxs^2) -lambda_s/2/dxs^2 ];
ATs_full(Ns,Nal+Np+Ns+2:Nal+Np+Ns+4)= [-lambda_s/2/dxs^2 (rho_s*Cps/dt+lambda_s/dxs^2) -lambda_s/2/dxs^2 ];

iTn = [2:Nn+1 2:Nn+1 2:Nn+1];
jTn = [1:Nn 2:Nn+1 3:Nn+2];
sTn = [-lambda_n/2/dxn^2*ones(1,Nn) (rho_n*Cpn/dt+lambda_n/dxn^2)*ones(1,Nn) -lambda_n/2/dxn^2*ones(1,Nn)];
ATn = sparse(iTn,jTn,sTn,Nn+2,Nn+2);ics = [2:Ns-1 2:Ns-1 2:Ns-1];
jcs = [1:Ns-2 2:Ns-1 3:Ns];
ATn_full = [sparse(Nn+2,Nal+Np+Ns+3) ATn sparse(Nn+2,Nco+1)];
ATn_full(1,Nal+Np+Ns+2:Nal+Np+Ns+6) = [lambda_s/dxs -4*lambda_s/dxs 3*(lambda_s/dxs+lambda_n/dxn) -4*lambda_n/dxn lambda_n/dxn];
ATn_full(Nn+2,Nal+Np+Ns+Nn+3:Nal+Np+Ns+Nn+7) = [lambda_n/dxn -4*lambda_n/dxn 3*(lambda_n/dxn+lambda_co/dxco) -4*lambda_co/dxco lambda_co/dxco];

iTco = [2:Nco 2:Nco 2:Nco];
jTco = [1:Nco-1 2:Nco 3:Nco+1];
sTco = [-lambda_co/2/dxco^2*ones(1,Nco-1) (rho_co*Cpco/dt+lambda_co/dxco^2)*ones(1,Nco-1) -lambda_co/2/dxco^2*ones(1,Nco-1)];
ATco = sparse(iTco,jTco,sTco,Nco+1,Nco+1);
ATco_full = [sparse(Nco+1,Nal+Np+Ns+Nn+5) ATco];
ATco_full(1,Nal+Np+Ns+Nn+5:Nal+Np+Ns+Nn+7)  = [-lambda_co/2/dxco^2 (rho_co*Cpco/dt+lambda_co/dxco^2) -lambda_co/2/dxco^2];
ATco_full(Nco+1,Nal+Np+Ns+Nn+Nco+4:Nal+Np+Ns+Nn+Nco+6) = [1 -4 3];

AT  = [ATal_full;ATp_full;ATs_full;ATn_full;ATco_full];
iAT = inv(AT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%From Initialize

icp = [2:Np+1 2:Np+1 2:Np+1];
jcp = [1:Np 2:Np+1 3:Np+2];

ics = [2:Ns-1 2:Ns-1 2:Ns-1];
jcs = [1:Ns-2 2:Ns-1 3:Ns];

icn = [2:Nn+1 2:Nn+1 2:Nn+1];
jcn = [1:Nn 2:Nn+1 3:Nn+2];


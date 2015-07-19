function f = T1potential(x)

global Prob1

%global Prob1.K.Nal Prob1.K.Np Prob1.K.Ns Prob1.K.Nn Prob1.K.Tref Prob1.K.Ip Prob1.K.In Prob1.K.Ini Prob1.K.Isi Prob1.K.Ipi Prob1.K.dp Prob1.K.dn Prob1.K.dxp Prob1.K.dxs Prob1.K.dxn Prob1.K.sig_eff Prob1.K.F Prob1.K.R Prob1.K.cs_max Prob1.K.i Prob1.K.dt Prob1.K.a_i Prob1.K.tplus  Prob1.K.Rp Prob1.M.iAc Prob1.M.iAp2 Prob1.M.iApphi1 Prob1.M.iAnphi1 

%Nal = Prob1.K.Prob1.K.Nal;
%Np  = Prob1.K.Prob1.K.Np;
%Ns  = Prob1.K.Prob1.K.Ns;
%Nn  = Prob1.K.Prob1.K.Nn;
%Tref= Prob1.K.Prob1.K.Tref;
%Ip  = Prob1.K.Prob1.K.Ip;
%In  = Prob1.K.Prob1.K.In;Prob1.M.I
%Ipi = Prob1.K.Prob1.K.Ipi;
%Isi = Prob1.K.Prob1.K.Isi;
%Ini = Prob1.K.Prob1.K.Ini;

%dp = Prob1.K.Prob1.K.dp;
%dn = Prob1.K.Prob1.K.dn;
%dxp= Prob1.K.Prob1.K.dxp;
%dxs=Prob1.K.Prob1.K.dxs;
%dxn=Prob1.K.Prob1.K.dxn;
%sig_eff = Prob1.K.Prob1.K.sig_eff;
%F  = Prob1.K.Prob1.K.F;
%R  = Prob1.K.Prob1.K.R;
%cs_max = Prob1.K.Prob1.K.cs_max;
%i = Prob1.K.Prob1.K.i;
%dt = Prob1.K.Prob1.K.dt;
%a_i = Prob1.K.Prob1.K.a_i;
%tplus = Prob1.K.Prob1.K.tplus;
%Rp = Prob1.K.Prob1.K.Rp;
%iAc = Prob1.M.Prob1.M.iAc;
%iAp2 = Prob1.M.Prob1.M.iAp2;
%iApphi1 = Prob1.M.Prob1.M.iApphi1;
%iAnphi1 = Prob1.M.Prob1.M.iAnphi1;
%k_pT = Prob1.K.Prob1.K.k_pT;
%k_nT = Prob1.K.Prob1.K.k_nT;

%jtemp = 1e-5*idct([[x(3:Prob1.K.dp+2);zeros(Prob1.K.Np+2-Prob1.K.dp,1)] [x(Prob1.K.dp+3:Prob1.K.dp+Prob1.K.dn+2);zeros(Prob1.K.Nn+2-Prob1.K.dn,1)]]);
%jtemp = 1e-5*idct([x(3:Prob1.K.dp+2) x(Prob1.K.dp+3:Prob1.K.dp+Prob1.K.dn+2)],Prob1.K.Np+2);
Prob1.M.j_i(Prob1.K.i+1,1:Prob1.K.Np+2) = 1e-5*Prob1.M.dpM*[x(3:Prob1.K.dp+2);zeros(Prob1.K.Np+2-Prob1.K.dp,1)];
Prob1.M.j_i(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4) = 1e-5*Prob1.M.dnM*[x(Prob1.K.dp+3:Prob1.K.dp+Prob1.K.dn+2);zeros(Prob1.K.Nn+2-Prob1.K.dn,1)];
% jtemp1 = 1e-5*spline(1:Prob1.K.dj:Prob1.K.Np+2,x(3:ceil((Prob1.K.Np+2)/Prob1.K.dj)+2),1:Prob1.K.Np+2);
% jtemp2 = 1e-5*spline(1:Prob1.K.dj:Prob1.K.Nn+2,x(ceil((Prob1.K.Np+2)/Prob1.K.dj)+3:end),1:Prob1.K.Nn+2);
%jtemp = 1e-5*x(3:end);
%Prob1.M.j_i(Prob1.K.i+1,[1:Prob1.K.Np+2 Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4]) = jtemp(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc  = [Prob1.M.fcp+[0 Prob1.K.a_i(1)*(1-Prob1.K.tplus(Prob1.K.Ip)).*Prob1.M.j_i(Prob1.K.i+1,Prob1.K.Ip) 0] Prob1.M.fcs Prob1.M.fcn+[0 Prob1.K.a_i(3)*(1-Prob1.K.tplus(Prob1.K.In)).*Prob1.M.j_i(Prob1.K.i+1,Prob1.K.In) 0]]';
Prob1.M.c(Prob1.K.i+1,:) = Prob1.M.iAc*fc;
Prob1.M.cs_avg(Prob1.K.i+1,:) = Prob1.M.cs_avg(Prob1.K.i,:) - 3*Prob1.K.dt/Prob1.K.Rp*Prob1.M.j_i(Prob1.K.i+1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpi = [0 Prob1.K.dxp*Prob1.K.a_i(1)*Prob1.K.F*Prob1.M.j_i(Prob1.K.i+1,Prob1.K.Ip) 2*Prob1.K.dxp*Prob1.K.a_i(1)*Prob1.K.F*Prob1.M.j_i(Prob1.K.i+1,Prob1.K.Np+2)]';
Prob1.M.fp1 = [x(1);Prob1.M.fp1(2:Prob1.K.Np+1);-2*Prob1.K.dxp*Prob1.M.I(Prob1.K.i+1)/Prob1.K.sig_eff(1)];

Prob1.M.Phi1_p = Prob1.M.iApphi1*[Prob1.M.fp1;fpi];

Prob1.M.Phi1(Prob1.K.i+1,1:Prob1.K.Np+2) = Prob1.M.Phi1_p(1:Prob1.K.Np+2)';
Prob1.M.ie(Prob1.K.i+1,1:Prob1.K.Np+2)   = Prob1.M.Phi1_p(Prob1.K.Np+3:end)';

fni = [Prob1.M.I(Prob1.K.i+1) Prob1.K.dxn*Prob1.K.a_i(3)*Prob1.K.F*Prob1.M.j_i(Prob1.K.i+1,Prob1.K.In) 2*Prob1.K.dxn*Prob1.K.a_i(3)*Prob1.K.F*Prob1.M.j_i(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4)]';
Prob1.M.fn1 = [x(2);Prob1.M.fn1(2:Prob1.K.Nn+1);-2*Prob1.K.dxn*Prob1.M.I(Prob1.K.i+1)/Prob1.K.sig_eff(3)];

Prob1.M.Phi1_n = Prob1.M.iAnphi1*[Prob1.M.fn1;fni];

Prob1.M.Phi1(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4) = Prob1.M.Phi1_n(1:Prob1.K.Nn+2)';
Prob1.M.ie(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4)   = Prob1.M.Phi1_n(Prob1.K.Nn+3:end)';

%%%%%%%%%%%%%%%%%%%%%%%%%
%Prob1.M.Phi2 reverse
%%%%%%%%%%%%%%%%%%%%%%%%%

Prob1.M.gp(2:Prob1.K.Nn+1,1)= ((Prob1.M.ie(Prob1.K.i+1,Prob1.K.Ini))*Prob1.K.dxp./(Prob1.V.Keff_n(Prob1.K.Nn:-1:1)) + Prob1.K.R*Prob1.M.T(Prob1.K.i,Prob1.K.Nal+Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4:-1:Prob1.K.Nal+Prob1.K.Np+Prob1.K.Ns+5)/Prob1.K.F.*(1-Prob1.K.tplus(Prob1.K.Ini))./(Prob1.M.c(Prob1.K.i+1,Prob1.K.Ini)).*(Prob1.M.c(Prob1.K.i+1,Prob1.K.Ini+1)-Prob1.M.c(Prob1.K.i+1,Prob1.K.Ini-1)));
Prob1.M.gp(Prob1.K.Nn+3:Prob1.K.Nn+Prob1.K.Ns+2,1)= ((Prob1.M.ie(Prob1.K.i+1,Prob1.K.Isi))*Prob1.K.dxs./(Prob1.V.Keff_s(Prob1.K.Ns:-1:1)) + Prob1.K.R*Prob1.M.T(Prob1.K.i,Prob1.K.Nal+Prob1.K.Np+Prob1.K.Ns+3:-1:Prob1.K.Nal+Prob1.K.Np+4)/Prob1.K.F.*(1-Prob1.K.tplus(Prob1.K.Isi))./(Prob1.M.c(Prob1.K.i+1,Prob1.K.Isi)).*(Prob1.M.c(Prob1.K.i+1,Prob1.K.Isi+1)-Prob1.M.c(Prob1.K.i+1,Prob1.K.Isi-1)));
Prob1.M.gp(Prob1.K.Nn+Prob1.K.Ns+4:Prob1.K.Nn+Prob1.K.Ns+Prob1.K.Np+3,1)= ((Prob1.M.ie(Prob1.K.i+1,Prob1.K.Ipi))*Prob1.K.dxn./(Prob1.V.Keff_p(Prob1.K.Np:-1:1)) + Prob1.K.R*Prob1.M.T(Prob1.K.i,Prob1.K.Nal+Prob1.K.Np+2:-1:Prob1.K.Nal+3)/Prob1.K.F.*(1-Prob1.K.tplus(Prob1.K.Ipi))./(Prob1.M.c(Prob1.K.i+1,Prob1.K.Ipi)).*(Prob1.M.c(Prob1.K.i+1,Prob1.K.Ipi+1)-Prob1.M.c(Prob1.K.i+1,Prob1.K.Ipi-1)));

Prob1.M.Phi2(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4:-1:1) = (Prob1.M.iAp2*Prob1.M.gp)';

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%
cs_surf(Prob1.K.i+1,1:Prob1.K.Np+2) = -Prob1.M.j_i(Prob1.K.i+1,1:Prob1.K.Np+2)/5*Prob1.K.Rp./Prob1.V.Dps_eff+Prob1.M.cs_avg(Prob1.K.i+1,1:Prob1.K.Np+2);
cs_surf(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4) = -Prob1.M.j_i(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4)/5*Prob1.K.Rp./Prob1.V.Dns_eff+Prob1.M.cs_avg(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4);

theta_p = cs_surf(Prob1.K.i+1,1:Prob1.K.Np+2)/Prob1.K.cs_max(1);
theta_n = cs_surf(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4)/Prob1.K.cs_max(3);

dudt_p = -0.001*(0.199521039-0.928373822*theta_p+1.364550689000003*theta_p.^2-0.6115448939999998*theta_p.^3);
dudt_p = dudt_p./(1-5.661479886999997*theta_p +11.47636191*theta_p.^2-9.82431213599998*theta_p.^3+3.048755063*theta_p.^4);

dudt_n = 0.001*(0.005269056 +3.299265709*theta_n-91.79325798*theta_n.^2+1004.911008*theta_n.^3-5812.278127*theta_n.^4 + ...
         19329.7549*theta_n.^5 - 37147.8947*theta_n.^6 + 38379.18127*theta_n.^7-16515.05308*theta_n.^8); % There is a typo in this eqn. in the original paper
dudt_n = dudt_n./(1-48.09287227*theta_n+1017.234804*theta_n.^2-10481.80419*theta_n.^3+59431.3*theta_n.^4-195881.6488*theta_n.^5+...
         374577.3152*theta_n.^6 - 385821.1607*theta_n.^7 + 165705.8597*theta_n.^8);
     
U_p   = (-4.656+88.669*theta_p.^2 - 401.119*theta_p.^4 + 342.909*theta_p.^6 - 462.471*theta_p.^8 + 433.434*theta_p.^10);
U_p   = U_p./(-1+18.933*theta_p.^2-79.532*theta_p.^4+37.311*theta_p.^6-73.083*theta_p.^8+95.96*theta_p.^10);
U_n   = 0.7222 + 0.1387*theta_n + 0.029*theta_n.^0.5 - 0.0172./theta_n + 0.0019./theta_n.^1.5 + 0.2808*exp(0.9-15*theta_n)-0.7984*exp(0.4465*theta_n - 0.4108);

U_p   = U_p + (Prob1.M.T(Prob1.K.i,Prob1.K.Nal+2:Prob1.K.Nal+Prob1.K.Np+3)-Prob1.K.Tref).*dudt_p;
U_n   = U_n + (Prob1.M.T(Prob1.K.i,Prob1.K.Nal+Prob1.K.Np+Prob1.K.Ns+4:Prob1.K.Nal+Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+5)-Prob1.K.Tref).*dudt_n;

deltap = (0.5*Prob1.K.F/Prob1.K.R./Prob1.M.T(Prob1.K.i,Prob1.K.Nal+2:Prob1.K.Nal+Prob1.K.Np+3).*(Prob1.M.Phi1(Prob1.K.i+1,1:Prob1.K.Np+2) - Prob1.M.Phi2(Prob1.K.i+1,1:Prob1.K.Np+2) - U_p));
deltan = (0.5*Prob1.K.F/Prob1.K.R./Prob1.M.T(Prob1.K.i,Prob1.K.Nal+Prob1.K.Np+Prob1.K.Ns+4:Prob1.K.Nal+Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+5).*(Prob1.M.Phi1(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4) - Prob1.M.Phi2(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4) - U_n));

kpinterim(1:Prob1.K.Np+2) = 2*Prob1.K.k_pT.*(Prob1.M.c(Prob1.K.i+1,1:Prob1.K.Np+2).*cs_surf(Prob1.K.i+1,1:Prob1.K.Np+2).*(Prob1.K.cs_max(1)-cs_surf(Prob1.K.i+1,1:Prob1.K.Np+2))).^0.5;
kpinterim(Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4) = 2*Prob1.K.k_nT.*(Prob1.M.c(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4).*cs_surf(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4).*(Prob1.K.cs_max(3)-cs_surf(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4))).^0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_1 = 10;
alpha_2 = 100;
alpha_3 = 100;
% alpha_4 = 0;
% alpha_5 = 0;

% if max(Prob1.M.ie(Prob1.K.i+1,:)) > 1
%     lambda_ie = 1000*exp(10*max(Prob1.M.ie(Prob1.K.i+1,:))-1);
% elseif min(Prob1.M.ie(Prob1.K.i+1,:)) < Prob1.M.I(Prob1.K.i+1,1)-1
%     lambda_ie = 1000*exp(10*min(-Prob1.M.ie(Prob1.K.i+1,:)) + Prob1.M.I(Prob1.K.i+1,1)+1);
% else
%     lambda_ie = 0;  
% end
% 
% if min(Prob1.M.c(Prob1.K.i+1,:)) < 0
%     lambda_ie1 = 100000*exp(-min(Prob1.M.c(Prob1.K.i+1,:)));
% else
%     lambda_ie1 = 0;
% end

% if max(deltap) > 0
%     lambda_ie2 = 0*exp(max(deltap));
% else
%     lambda_ie2 = 0;
% end
% 
% if min(deltan) < 0
%     lambda_ie3 = 0*exp(max(deltan));
% else
%     lambda_ie3 = 0;
% end

f = alpha_1*(norm(([Prob1.M.ie(Prob1.K.i+1,Prob1.K.Np+2);Prob1.M.ie(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4)]-[Prob1.M.I(Prob1.K.i+1);0]))) + ...
    alpha_2*norm(deltap-asinh(Prob1.M.j_i(Prob1.K.i+1,1:Prob1.K.Np+2)./kpinterim(1:Prob1.K.Np+2))) + ...
    alpha_3*norm(deltan-asinh(Prob1.M.j_i(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4)./kpinterim(Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4)));%+...
%     alpha_4*norm(1e5*diff(Prob1.M.j_i(Prob1.K.i+1,1:Prob1.K.Np+2)))+...
%     alpha_5*norm(1e5*diff(Prob1.M.j_i(Prob1.K.i+1,Prob1.K.Np+Prob1.K.Ns+3:Prob1.K.Np+Prob1.K.Ns+Prob1.K.Nn+4)));

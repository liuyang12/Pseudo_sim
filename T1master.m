%This program models the psuedo 2D Lithium ion battery model from the
%artice :
%Simple finite difference is used to numerically solve the PDEs.
clear all
tic

global Prob1

Time = 5000; % Duration of discharge/charge
Nt= 500; % Number of time intervals
% [al] aluminum current collector
% [p]  positive electrode
% [s]  separator
% [n]  negative electrode
% [co] copper current collector

Nal=5;      % 
Np= 25;     % Number of spatial intervals on the positive electrode side
Ns= 5;      % Number of spatial intervals in the separator
Nn= 25;     % Number of spatial intervals on the negative electrode side
Nco=5;      % 

dt = Time/Nt; % Size of time interval in secs

%--------------------------------------------------------------------------

T1params; % Program with constants

%--------------------------------------------------------------------------
% MEMORY Initialization

fTal= zeros(1,Nal+1);
Prob1.M.fcp = zeros(1,Np+2);
Prob1.M.fcn = zeros(1,Nn+2);
fTco= zeros(1,Nco+1);
phi0_vec = zeros(Nt,2);
V = zeros(Nt,2);

Prob1.M.c    = zeros(Nt,Np+Ns+Nn+4);
Prob1.M.Phi1 = zeros(Nt,Np+Ns+Nn+4);
Prob1.M.Phi2 = zeros(Nt,Np+Ns+Nn+4);
Prob1.M.j_i  = zeros(Nt,Np+Ns+Nn+4);
Prob1.M.ie   = zeros(Nt,Np+Ns+Nn+4);
Prob1.M.T    = zeros(Nt,Nal+Np+Ns+Nn+Nco+6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prob1.M.I = idinput(Nt+2,'rbs',[0 0.1],[-50 -20]);%-30;
Prob1.M.I = -30*ones(Time,1);
TempUpdate = 1;
phip0= 4.156;
phin0= 0.0747;
V(1:2)     = 4.2;
disp('Discharging ...  % 00');
bb    = 3;

i = 0;
dp = 2;
dn = 2;
Prob1.K.dj = 3;
Prob1.K.Nal = Nal;
Prob1.K.Np  = Np;
Prob1.K.Ns  = Ns;
Prob1.K.Nn  = Nn;
Prob1.K.Tref= Tref;
Prob1.K.Ip  = Ip;
Prob1.K.In  = In;
Prob1.K.Ipi = Ipi;
Prob1.K.Isi = Isi;
Prob1.K.Ini = Ini;

Prob1.K.dp  = dp;
Prob1.K.dn  = dn;
Prob1.K.dxp = dxp;
Prob1.K.dxs = dxs;
Prob1.K.dxn = dxn;
Prob1.K.sig_eff = sig_eff;
Prob1.K.F = F;
Prob1.K.R = R;
Prob1.K.cs_max = cs_max;
Prob1.K.i = i;
Prob1.K.dt = dt;
Prob1.K.a_i = a_i;
Prob1.K.tplus = tplus*ones(1,Np+Ns+Nn+4);
Prob1.K.Rp = Rp;
Prob1.M.Phi1 = [phip0*ones(1,Np+2) zeros(1,Ns) phin0*ones(1,Nn+2)];
Prob1.M.dpM = inv(dctmtx(Np+2));
Prob1.M.dnM = inv(dctmtx(Nn+2));

while (V(i+1) >= 2.6) && (i <= Nt)
    i = i+1;
    s = sprintf('\b');
    s = repmat(s,1,bb);
    fprintf(s);
    s = sprintf('%d%%',ceil(100/Nt*i));
    fprintf(s);
    bb = length(s)-1;
    
    T1initialize; % Concentration Calculations
    
    if i == 1
        
        cp = 1e5*dct([Prob1.M.j_i(i+1,1:Np+2)' Prob1.M.j_i(i+1,Np+Ns+3:Np+Ns+Nn+4)']);
        LB = [-inf;-inf;-inf*ones(dp+dn,1)];%[-inf;-inf;-inf*ones(Np+2,1);zeros(Nn+2,1)];
        LU = [inf;inf;inf*ones(dp+dn,1)];%[inf;inf;zeros(Np+2,1);inf*ones(Nn+2,1)];
%         cp = 1e5*([Prob1.M.j_i(i+1,1:Np+2)' Prob1.M.j_i(i+1,Np+Ns+3:Np+Ns+Nn+4)']);
%         LB = [-inf;-inf;-inf*ones(Np+2,1);zeros(Nn+2,1)];
%         LU = [inf;inf;zeros(Np+2,1);inf*ones(Nn+2,1)];
        x_opt = [];
        f_opt = [];
        f_Low = [];
        x_min = [];
        x_max = [];
        Name = 'test';
        x_0= [phip0;phin0;cp(1:dp,1);cp(1:dn,2)]';
        
        Prob.KNITRO.options.ALG = 3;
        Acon = [eye(2) zeros(2,dp+dn);zeros(Np+2,2) Prob1.M.dpM(:,1:dp) zeros(Np+2,dn);zeros(Nn+2,2+dp) Prob1.M.dnM(:,1:dn)];%*[x(3:Prob1.K.dp+2)];
        b_L = [0*ones(2,1);-inf*ones(Np+2,1);zeros(Nn+2,1)];
        b_U = [10*ones(2,1);zeros(Np+2,1);inf*ones(Nn+2,1)];

        Prob = conAssign('T1potential',[],[],[],LB,LU,Name,x_0,[],0,Acon,b_L,b_U);
        Prob.KNITRO.options.ALG = 3;
        %Prob.A = [ones(2,2+dp+dn);zeros(Np+2,2) Prob1.M.dpM(:,1:dp) zeros(Np+2,dn);zeros(Nn+2,2+dp) Prob1.M.dnM(:,1:dn)];%*[x(3:Prob1.K.dp+2)];
%         Prob.A = eye(2+dp+dn);
%         
%         Prob.b_L = [zeros(2,1);-inf*ones(dp,1);zeros(dn,1)];
%         Prob.b_U = [zeros(2,1);zeros(dp,1);inf*ones(dn,1)];
    else
        
        x_0= [phip0;phin0;cp(1:dp,1);cp(1:dn,2)]';
        
        Acon = [eye(2) zeros(2,dp+dn);zeros(Np+2,2) Prob1.M.dpM(:,1:dp) zeros(Np+2,dn);zeros(Nn+2,2+dp) Prob1.M.dnM(:,1:dn)];%*[x(3:Prob1.K.dp+2)];
        b_L = [-10*ones(2,1);-inf*ones(Np+2,1);zeros(Nn+2,1)];
        b_U = [10*ones(2,1);zeros(Np+2,1);inf*ones(Nn+2,1)];

        Prob = conAssign('T1potential',[],[],[],LB,LU,Name,x.x_k,[],0,Acon,b_L,b_U);
        Prob.KNITRO.options.ALG = 3;
        %Prob.A = [ones(2,2+dp+dn);zeros(Np+2,2) Prob1.M.dpM(:,1:dp) zeros(Np+2,dn);zeros(Nn+2,2+dp) Prob1.M.dnM(:,1:dn)];%*[x(3:Prob1.K.dp+2)];
 
        %Prob.b_L = [zeros(2,1);-inf*ones(Np+2,1);zeros(Nn+2,1)];
        %Prob.b_U = [zeros(2,1);zeros(Np+2,1);inf*ones(Nn+2,1)];
     end
    Prob.KNITRO.options.FEASTOL = 0.1;
    Prob.KNITRO.options.FEASTOL_ABS = 0.1;
    Prob.KNITRO.options.MAXIT = 10;
    Prob.KNITRO.options.OPTTOL = 0.1;
    Prob.KNITRO.options.OPTTOL_ABS = 0.1;
    Prob1.K.i = i;

    x = tomRun('knitro',Prob,0);
    Prob      = WarmDefSOL('knitro',Prob,x);
    phip0 = x.x_k(1);
    phin0 = x.x_k(2);
    jtemp = 1e-5*idct([[x.x_k(3:dp+2);zeros(Np+2-dp,1)] [x.x_k(dp+3:dp+dn+2);zeros(Nn+2-dn,1)]]);
 %   jtemp = 1e-5*x.x_k(3:end);
%     jtemp1 = 1e-5*spline(1:Prob1.K.dj:Prob1.K.Np+2,x.x_k(3:ceil((Prob1.K.Np+2)/Prob1.K.dj)+2),1:Prob1.K.Np+2);
%     jtemp2 = 1e-5*spline(1:Prob1.K.dj:Prob1.K.Nn+2,x.x_k(ceil((Prob1.K.Np+2)/Prob1.K.dj)+3:end),1:Prob1.K.Nn+2);
% 
%     
    Prob1.M.j_i(i+1,[1:Np+2 Np+Ns+3:Np+Ns+Nn+4]) = jtemp(:)';
    
    T1Temp;
    if i >=2
        
        Prob1.M.fcp(2:Np+1) = Deff_p(i+1,2:Np+1)/2/dxp^2.*(Prob1.M.c(i,Ip+1)+Prob1.M.c(i,Ip-1)) + (eps_p/dt-Deff_p(i+1,2:Np+1)/dxp^2).*Prob1.M.c(i,Ip) + a_i(1)*(1-tplus)*Prob1.M.j_i(i+1,Ip);
        
        Prob1.M.fcs = Deff_s(i+1,2:end-1)/2/dxs^2.*(Prob1.M.c(i,Is+1)+Prob1.M.c(i,Is-1)) + (eps_s/dt-Deff_s(i+1,2:end-1)/dxs^2).*Prob1.M.c(i,Is) + a_i(2)*(1-tplus)*Prob1.M.j_i(i+1,Is);
        
        Prob1.M.fcn(2:Nn+1) = Deff_n(i+1,2:Nn+1)/2/dxn^2.*(Prob1.M.c(i,In+1)+Prob1.M.c(i,In-1)) + (eps_n/dt-Deff_n(i+1,2:Nn+1)/dxn^2).*Prob1.M.c(i,In) + a_i(3)*(1-tplus)*Prob1.M.j_i(i+1,In);
        fc  = [Prob1.M.fcp';Prob1.M.fcs';Prob1.M.fcn'];
        
        Prob1.M.c(i+1,:) = Ac\fc;
        Prob1.M.cs_avg(i+1,:) = Prob1.M.cs_avg(i,:) - 3*dt/Rp*Prob1.M.j_i(i+1,:);
        
    end
    if mod(i,5) == 0
        %keyboard
    end
    phi0_vec(i,:) =[phip0 phin0];
    
%save T_master    
end
etime = toc;
T = Prob1.M.T;
c = Prob1.M.c;
ie= Prob1.M.ie;
Phi1=Prob1.M.Phi1;
Phi2=Prob1.M.Phi2;
j_i=Prob1.M.j_i;

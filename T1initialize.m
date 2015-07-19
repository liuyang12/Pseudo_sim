
if i == 1
    
    Prob1.M.c(1,1:Np+Ns+Nn+4) = 1000*ones(1,Np+Ns+Nn+4);
    Prob1.M.c(2,1:Np+Ns+Nn+4) = 1000*ones(1,Np+Ns+Nn+4);
    
    Prob1.M.cs_avg(1,1:Np+2) = cs_initp*ones(1,Np+2);
    Prob1.M.cs_avg(1,Np+3:Np+Ns+2) = cs_inits*ones(1,Ns);
    Prob1.M.cs_avg(1,Np+Ns+3:Np+Ns+Nn+4) = cs_initn*ones(1,Nn+2);
    
    Prob1.M.cs_avg(2,1:Np+2) = cs_initp*ones(1,Np+2);
    Prob1.M.cs_avg(2,Np+3:Np+Ns+2) = cs_inits*ones(1,Ns);
    Prob1.M.cs_avg(2,Np+Ns+3:Np+Ns+Nn+4) = cs_initn*ones(1,Nn+2);
    
    Prob1.M.T(1,1:Nal+Np+Ns+Nn+Nco+6) = Tref;
    Prob1.M.T(2,1:Nal+Np+Ns+Nn+Nco+6) = Tref;
   
else

    if i >= 1
        
        V(i+1) = [Prob1.M.Phi1(i,1)-Prob1.M.Phi1(i,end)];
        
    end
    
end

if i == 1
    Prob1.M.j_i(i+1,:) = [-4.3916e-6*ones(1,Np+2) zeros(1,Ns) 4.8826e-6*ones(1,Nn+2)];
    Prob1.M.ie(i+1,:)  = [(0:Prob1.M.I(i+1)/Np:Prob1.M.I(i+1)) Prob1.M.I(i+1) Prob1.M.I(i+1)*ones(1,Ns) Prob1.M.I(i+1):-Prob1.M.I(i+1)/Nn:0 0];
else
    Prob1.M.j_i(i+1,:) = Prob1.M.j_i(i,:);
    Prob1.M.ie(i+1,:)  = [(0:Prob1.M.I(i+1)/Np:Prob1.M.I(i+1)) Prob1.M.I(i+1) Prob1.M.I(i+1)*ones(1,Ns) Prob1.M.I(i+1):-Prob1.M.I(i+1)/Nn:0 0];
end
%%%%%%%%%%%%%%
% Prob1.K.tplus = [ones(length(Prob1.M.c(i,:)),1) Prob1.M.c(i,:).^(0.5)' Prob1.M.c(i,:).^(1.5)' Prob1.M.T(i,Nal+2:Nal+Np+Ns+Nn+5).^(0.5)'.*Prob1.M.c(i,:).^(1.5)' Prob1.M.T(i,Nal+2:Nal+Np+Ns+Nn+5)'.*Prob1.M.c(i,:).^(1.5)' Prob1.M.T(i,Nal+2:Nal+Np+Ns+Nn+5).^(1.5)'.*Prob1.M.c(i,:).^(1.5)'];
% temptplus = [1.0133;-0.0128;0.5280;-0.0904;0.0052;-0.0001];
% Prob1.K.tplus = [Prob1.K.tplus*temptplus]';
%%%%%%%%%%%%%%%

Prob1.M.fp1 = -dxp*Prob1.M.I(i+1)/sig_eff(1)*ones(Np+2,1);
Prob1.M.Phi1(i+1,1:Np+Ns+Nn+4) = zeros(1,Np+Ns+Nn+4);

Prob1.M.ie(i+1,Np+3:Np+Ns+2) = Prob1.M.I(i+1)*ones(1,Ns);
Prob1.M.fn1 = -dxn*Prob1.M.I(i+1)/sig_eff(3)*ones(Nn+2,1);

%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%Diffusion Coefficients
%%%%%%%%%%%%%%%%%
% if TempUpdate == 0
%     if i == 1
%         Deff_p(i+1,:) = eps_p^brugg_p*1e-4*10.^((-4.43-54./(Prob1.M.T(i,Nal+2:Nal+Np+3)-229-5e-3*Prob1.M.c(i,1:Np+2))-0.22e-3*Prob1.M.c(i,1:Np+2)));
%         Deff_s(i+1,:) = eps_s^brugg_s*1e-4*10.^((-4.43-54./(Prob1.M.T(i,Nal+Np+3:Nal+Np+Ns+4)-229-5e-3*Prob1.M.c(i,Np+2:Np+Ns+3))-0.22e-3*Prob1.M.c(i,Np+2:Np+Ns+3)));
%         Deff_n(i+1,:) = eps_n^brugg_n*1e-4*10.^((-4.43-54./(Prob1.M.T(i,Nal+Np+Ns+4:Nal+Np+Ns+Nn+5)-229-5e-3*Prob1.M.c(i,Np+Ns+3:Np+Ns+Nn+4))-0.22e-3*Prob1.M.c(i,Np+Ns+3:Np+Ns+Nn+4)));
%         
%         Prob1.V.Dps_eff= Dps*exp(-EaDs/R*(1./Prob1.M.T(i,Nal+2:Nal+Np+3)-1/Tref));
%         Prob1.V.Dns_eff= Dns*exp(-Eaki/R*(1./Prob1.M.T(i,Nal+Np+Ns+4:Nal+Np+Ns+Nn+5)-1/Tref));
%         
%         Prob1.K.k_pT = k_p*exp(-5000/R*(1./Prob1.M.T(i,Nal+2:Nal+Np+3)-1/Tref));
%         Prob1.K.k_nT = k_n*exp(-5000/R*(1./Prob1.M.T(i,Nal+Np+Ns+4:Nal+Np+Ns+Nn+5)-1/Tref));
%         
%         Keff(i+1,:) = 1e-4*Prob1.M.c(i,:).*((-10.5+0.668*1e-3*Prob1.M.c(i,:)+0.494*1e-6*Prob1.M.c(i,:).^2) +...
%             (0.074 -0.0178*1e-3*Prob1.M.c(i,:)-8.86e-4*1e-6*Prob1.M.c(i,:).^2).*Prob1.M.T(i,Nal+2:end-Nco-1) + (-6.96e-5+2.8e-5*1e-3*Prob1.M.c(i,:)).*Prob1.M.T(i,Nal+2:end-Nco-1).^2).^2;
%     else
%         Deff_p(i+1,:) = Deff_p(i,:);
%         Deff_s(i+1,:) = Deff_s(i,:);
%         Deff_n(i+1,:) = Deff_n(i,:);
%         Keff(i+1,:)   = Keff(i,:);
%     end
% else
    Deff_p(i+1,:) = eps_p^brugg_p*1e-4*10.^((-4.43-54./(Prob1.M.T(i,Nal+2:Nal+Np+3)-229-5e-3*Prob1.M.c(i,1:Np+2))-0.22e-3*Prob1.M.c(i,1:Np+2)));
    Deff_s(i+1,:) = eps_s^brugg_s*1e-4*10.^((-4.43-54./(Prob1.M.T(i,Nal+Np+3:Nal+Np+Ns+4)-229-5e-3*Prob1.M.c(i,Np+2:Np+Ns+3))-0.22e-3*Prob1.M.c(i,Np+2:Np+Ns+3)));
    Deff_n(i+1,:) = eps_n^brugg_n*1e-4*10.^((-4.43-54./(Prob1.M.T(i,Nal+Np+Ns+4:Nal+Np+Ns+Nn+5)-229-5e-3*Prob1.M.c(i,Np+Ns+3:Np+Ns+Nn+4))-0.22e-3*Prob1.M.c(i,Np+Ns+3:Np+Ns+Nn+4)));
    
    Prob1.V.Dps_eff= Dps*exp(-EaDs/R*(1./Prob1.M.T(i,Nal+2:Nal+Np+3)-1/Tref));
    Prob1.V.Dns_eff= Dns*exp(-EaDs/R*(1./Prob1.M.T(i,Nal+Np+Ns+4:Nal+Np+Ns+Nn+5)-1/Tref));
    
    Prob1.K.k_pT = k_p*exp(-Eaki/R*(1./Prob1.M.T(i,Nal+2:Nal+Np+3)-1/Tref));
    Prob1.K.k_nT = k_n*exp(-Eaki/R*(1./Prob1.M.T(i,Nal+Np+Ns+4:Nal+Np+Ns+Nn+5)-1/Tref));
    
    Keff(i+1,:) = 1e-4/Prob1.K.C*Prob1.M.c(i,:).*((-10.5+0.668*1e-3*Prob1.M.c(i,:)+0.494*1e-6*Prob1.M.c(i,:).^2) +...
        (0.074 -0.0178*1e-3*Prob1.M.c(i,:)-8.86e-4*1e-6*Prob1.M.c(i,:).^2).*Prob1.M.T(i,Nal+2:end-Nco-1) + (-6.96e-5+2.8e-5*1e-3*Prob1.M.c(i,:)).*Prob1.M.T(i,Nal+2:end-Nco-1).^2).^2;
    
% end

Prob1.V.Keff_p = eps_p^brugg_p*Keff(i+1,Ip);
Prob1.V.Keff_s = eps_s^brugg_s*Keff(i+1,Is);
Prob1.V.Keff_n = eps_n^brugg_n*Keff(i+1,In);

Prob1.V.Keff_p0 = eps_p^brugg_p*Keff(i+1,1);
Prob1.V.Keff_p1 = eps_p^brugg_p*Keff(i+1,Np+2);
Prob1.V.Keff_s1 = eps_s^brugg_s*Keff(i+1,Np+2);
Prob1.V.Keff_s2 = eps_s^brugg_s*Keff(i+1,Np+Ns+3);
Prob1.V.Keff_n0 = eps_n^brugg_n*Keff(i+1,Np+Ns+3);
Prob1.V.Keff_n1 = eps_n^brugg_n*Keff(i+1,Np+Ns+Nn+4);
%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concentration
%%%%%%%%%%%%%%%%%%%%%%%%%%%

scp = [-Deff_p(i+1,2:Np+1)/2/dxp^2 (eps_p/dt+Deff_p(i+1,2:Np+1)/dxp^2) -Deff_p(i+1,2:Np+1)/2/dxp^2];
Acp = sparse(icp,jcp,scp,Np+2,Np+2);
Acp(1,1:3)  = [-3 4 -1];
Acp_full = [Acp sparse(Np+2,Ns+Nn+2)];
Acp_full(Np+2,Np:Np+4) = [Deff_p(i+1,Np)/dxp -4*Deff_p(i+1,Np+1)/dxp 3*(Deff_p(i+1,Np+2)/dxp+Deff_s(i+1,1)/dxs) -4*Deff_s(i+1,2)/dxs Deff_s(i+1,3)/dxs];

scs = [-Deff_s(i+1,3:end-2)/2/dxs^2 (eps_s/dt+Deff_s(i+1,3:end-2)/dxs^2) -Deff_s(i+1,3:end-2)/2/dxs^2];
Acs = sparse(ics,jcs,scs,Ns,Ns);
Acs_full = [sparse(Ns,Np+2) Acs sparse(Ns,Nn+2)];
Acs_full(1,Np+2:Np+4) = [-Deff_s(i+1,1)/2/dxs^2 (eps_s/dt+Deff_s(i+1,2)/dxs^2) -Deff_s(i+1,3)/2/dxs^2 ];
Acs_full(Ns,Np+Ns+1:Np+Ns+3)= [-Deff_s(i+1,end-2)/2/dxs^2 (eps_s/dt+Deff_s(i+1,end-1)/dxs^2) -Deff_s(i+1,end)/2/dxs^2 ];

scn = [-Deff_n(i+1,2:Nn+1)/2/dxn^2 (eps_n/dt+Deff_n(i+1,2:Nn+1)/dxn^2) -Deff_n(i+1,2:Nn+1)/2/dxn^2];
Acn = sparse(icn,jcn,scn,Nn+2,Nn+2);
Acn(end,end-2:end) = [1 -4 3];
Acn_full = [sparse(Nn+2,Np+Ns+2) Acn];
Acn_full(1,Np+Ns+1:Np+Ns+5) = [Deff_s(i+1,end-2)/dxs -4*Deff_s(i+1,end-1)/dxs 3*(Deff_s(i+1,end)/dxs+Deff_n(i+1,1)/dxn) -4*Deff_n(i+1,2)/dxn Deff_n(i+1,3)/dxn];
Ac  = [Acp_full;Acs_full;Acn_full];
Prob1.M.iAc = inv(Ac);

Prob1.M.fcp(2:Np+1) = Deff_p(i+1,2:Np+1)/2/dxp^2.*(Prob1.M.c(i,Ip+1)+Prob1.M.c(i,Ip-1)) + (eps_p/dt-Deff_p(i+1,2:Np+1)/dxp^2).*Prob1.M.c(i,Ip);
Prob1.M.fcs = Deff_s(i+1,2:end-1)/2/dxs^2.*(Prob1.M.c(i,Is+1)+Prob1.M.c(i,Is-1)) + (eps_s/dt-Deff_s(i+1,2:end-1)/dxs^2).*Prob1.M.c(i,Is);
Prob1.M.fcn(2:Nn+1) = Deff_n(i+1,2:Nn+1)/2/dxn^2.*(Prob1.M.c(i,In+1)+Prob1.M.c(i,In-1)) + (eps_n/dt-Deff_n(i+1,2:Nn+1)/dxn^2).*Prob1.M.c(i,In);
%%%%%%%%%%%%%%%%%%
%Prob1.M.Phi1
%%%%%%%%%%%%%%%%%%
Prob1.M.fp1 = [phip0;Prob1.M.fp1(2:Np+1);-2*dxp*Prob1.M.I(i+1)/sig_eff(1)];
Prob1.M.fn1 = [phin0;Prob1.M.fn1(2:Nn+1);-2*dxn*Prob1.M.I(i+1)/sig_eff(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ap2(Nn+2,Nn:Nn+4) = [Prob1.V.Keff_n1/dxn -4*Prob1.V.Keff_n1/dxn 3*(Prob1.V.Keff_n1/dxn+Prob1.V.Keff_s2/dxs) -4*Prob1.V.Keff_s2/dxs Prob1.V.Keff_s2/dxs];
Ap2(Nn+Ns+3,Nn+Ns+1:Nn+Ns+5) = [Prob1.V.Keff_s1/dxs -4*Prob1.V.Keff_s1/dxs 3*(Prob1.V.Keff_s1/dxs+Prob1.V.Keff_p1/dxp) -4*Prob1.V.Keff_p1/dxp Prob1.V.Keff_p1/dxp];
Prob1.M.iAp2 = inv(Ap2);

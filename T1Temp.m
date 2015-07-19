
% Keff(i+1,:) = 1e-4*Prob1.M.c(i+1,:).*((-10.5+0.668*1e-3*Prob1.M.c(i+1,:)+0.494*1e-6*Prob1.M.c(i+1,:).^2) +...
%     (0.074 -0.0178*1e-3*Prob1.M.c(i+1,:)-8.86e-4*1e-6*Prob1.M.c(i+1,:).^2).*Prob1.M.T(i,Nal+2:end-Nco-1) + (-6.96e-5+2.8e-5*1e-3*Prob1.M.c(i+1,:)).*Prob1.M.T(i,Nal+2:end-Nco-1).^2).^2;
% Keff_p = eps_p^brugg_p*Keff(i+1,Ip);
% Keff_s = eps_s^brugg_s*Keff(i+1,Is);
% Keff_n = eps_n^brugg_n*Keff(i+1,In);
%
% Keff_p0 = eps_p^brugg_p*Keff(i+1,1);
% Keff_p1 = eps_p^brugg_p*Keff(i+1,Np+2);
% Keff_s1 = eps_s^brugg_s*Keff(i+1,Np+2);
% Keff_s2 = eps_s^brugg_s*Keff(i+1,Np+Ns+3);
% Keff_n0 = eps_n^brugg_n*Keff(i+1,Np+Ns+3);
% Keff_n1 = eps_n^brugg_n*Keff(i+1,Np+Ns+Nn+4);

if TempUpdate == 0
    
    Prob1.M.T(i+1,:) = Tref;
    
else
    
    cs_surf(i+1,1:Np+2) = -Prob1.M.j_i(i+1,1:Np+2)/5*Rp./Prob1.V.Dps_eff+Prob1.M.cs_avg(i+1,1:Np+2);
    cs_surf(i+1,Np+Ns+3:Np+Ns+Nn+4) = -Prob1.M.j_i(i+1,Np+Ns+3:Np+Ns+Nn+4)/5*Rp./Prob1.V.Dns_eff+Prob1.M.cs_avg(i+1,Np+Ns+3:Np+Ns+Nn+4);
    
    theta_p = cs_surf(i+1,1:Np+2)/cs_max(1);
    theta_n = cs_surf(i+1,Np+Ns+3:Np+Ns+Nn+4)/cs_max(3);
    
    dudt_p = -0.001*(0.199521039-0.928373822*theta_p+1.364550689000003*theta_p.^2-0.6115448939999998*theta_p.^3);
    dudt_p = dudt_p./(1-5.661479886999997*theta_p +11.47636191*theta_p.^2-9.82431213599998*theta_p.^3+3.048755063*theta_p.^4);
    
    dudt_n = 0.001*(0.005269056 +3.299265709*theta_n-91.79325798*theta_n.^2+1004.911008*theta_n.^3-5812.278127*theta_n.^4 + ...
        19329.7549*theta_n.^5 - 37147.8947*theta_n.^6 + 38379.18127*theta_n.^7-16515.05308*theta_n.^8); % There is a typo in this eqn. in the original paper
    dudt_n = dudt_n./(1-48.09287227*theta_n+1017.234804*theta_n.^2-10481.80419*theta_n.^3+59431.3*theta_n.^4-195881.6488*theta_n.^5+...
        374577.3152*theta_n.^6 - 385821.1607*theta_n.^7 + 165705.8597*theta_n.^8);
    
    U_p   = (-4.656+88.669*theta_p.^2 - 401.119*theta_p.^4 + 342.909*theta_p.^6 - 462.471*theta_p.^8 + 433.434*theta_p.^10);
    U_p   = U_p./(-1+18.933*theta_p.^2-79.532*theta_p.^4+37.311*theta_p.^6-73.083*theta_p.^8+95.96*theta_p.^10);
    U_n   = 0.7222 + 0.1387*theta_n + 0.029*theta_n.^0.5 - 0.0172./theta_n + 0.0019./theta_n.^1.5 + 0.2808*exp(0.9-15*theta_n)-0.7984*exp(0.4465*theta_n - 0.4108);
    
    U_p   = U_p + (Prob1.M.T(i,Nal+2:Nal+Np+3)-Tref).*dudt_p;
    U_n   = U_n + (Prob1.M.T(i,Nal+Np+Ns+4:Nal+Np+Ns+Nn+5)-Tref).*dudt_n;
    
    deltap = (0.5*F/R./Prob1.M.T(i,Nal+2:Nal+Np+3).*(Prob1.M.Phi1(i+1,1:Np+2) - Prob1.M.Phi2(i+1,1:Np+2) - U_p));
    deltan = (0.5*F/R./Prob1.M.T(i,Nal+Np+Ns+4:Nal+Np+Ns+Nn+5).*(Prob1.M.Phi1(i+1,Np+Ns+3:Np+Ns+Nn+4) - Prob1.M.Phi2(i+1,Np+Ns+3:Np+Ns+Nn+4) - U_n));
    
    kpinterim(1:Np+2) = 2*Prob1.K.k_pT.*(Prob1.M.c(i+1,1:Np+2).*cs_surf(i+1,1:Np+2).*(cs_max(1)-cs_surf(i+1,1:Np+2))).^0.5;
    kpinterim(Np+Ns+3:Np+Ns+Nn+4) = 2*Prob1.K.k_nT.*(Prob1.M.c(i+1,Np+Ns+3:Np+Ns+Nn+4).*cs_surf(i+1,Np+Ns+3:Np+Ns+Nn+4).*(cs_max(3)-cs_surf(i+1,Np+Ns+3:Np+Ns+Nn+4))).^0.5;
    
    ji_int(1:Np+2) = kpinterim(1:Np+2).*sinh(deltap);
    ji_int(Np+Ns+3:Np+Ns+Nn+4) = kpinterim(Np+Ns+3:Np+Ns+Nn+4).*sinh(deltan);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Qrxn_p = F*a_i(1)*Prob1.M.j_i(i+1,1:Np+2).*(Prob1.M.Phi1(i+1,1:Np+2) - Prob1.M.Phi2(i+1,1:Np+2) - U_p);
    Qrxn_n = F*a_i(3)*Prob1.M.j_i(i+1,Np+Ns+3:Np+Ns+Nn+4).*(Prob1.M.Phi1(i+1,Np+Ns+3:Np+Ns+Nn+4) - Prob1.M.Phi2(i+1,Np+Ns+3:Np+Ns+Nn+4) - U_n);
    
    Qrev_p = F*a_i(1)*Prob1.M.j_i(i+1,Np+Ns+3:Np+Ns+Nn+4).*Prob1.M.T(i,Nal+2:Nal+Np+3).*dudt_p;
    Qrev_n = F*a_i(3)*Prob1.M.j_i(i+1,Np+Ns+3:Np+Ns+Nn+4).*Prob1.M.T(i,Nal+Np+Ns+4:Nal+Np+Ns+Nn+5).*dudt_n;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Qohm_p(1,Ip) = sig_eff(1)*((Prob1.M.Phi1(i+1,Ip+1)-Prob1.M.Phi1(i+1,Ip-1))/2/dxp).^2+Prob1.V.Keff_p.*((Prob1.M.Phi2(i+1,Ip+1)-Prob1.M.Phi2(i+1,Ip-1))/2/dxp).^2 + ...
        2*Prob1.V.Keff_p*R.*Prob1.M.T(i,Nal+3:Nal+Np+2)/F.*(1-Prob1.K.tplus(Ip))./Prob1.M.c(i+1,Ip).*(Prob1.M.c(i+1,Ip+1)-Prob1.M.c(i+1,Ip-1))/2/dxp.*(Prob1.M.Phi2(i+1,Ip+1)-Prob1.M.Phi2(i+1,Ip-1))/2/dxp;
    Qohm_n(1,2:Nn+1) = sig_eff(3)*((Prob1.M.Phi1(i+1,In+1)-Prob1.M.Phi1(i+1,In-1))/2/dxn).^2+Prob1.V.Keff_n.*((Prob1.M.Phi2(i+1,In+1)-Prob1.M.Phi2(i+1,In-1))/2/dxn).^2 + ...
        2*Prob1.V.Keff_n*R.*Prob1.M.T(i,Nal+Np+Ns+5:Nal+Np+Ns+Nn+4)/F.*(1-Prob1.K.tplus(In))./Prob1.M.c(i+1,In).*(Prob1.M.c(i+1,In+1)-Prob1.M.c(i+1,In-1))/2/dxn.*(Prob1.M.Phi2(i+1,In+1)-Prob1.M.Phi2(i+1,In-1))/2/dxn;
    
    Qohm_p(1)   = sig_eff(1)*((-3*Prob1.M.Phi1(i+1,1)+4*Prob1.M.Phi1(i+1,2)-Prob1.M.Phi1(i+1,3))/2/dxp).^2+Prob1.V.Keff_p0.*...
        ((-3*Prob1.M.Phi2(i+1,1)+4*Prob1.M.Phi2(i+1,2)-Prob1.M.Phi2(i+1,3))/2/dxp).^2 + ...
        2*Prob1.V.Keff_p0*R.*Prob1.M.T(i,Nal+2)/F.*(1-Prob1.K.tplus(1))./Prob1.M.c(i+1,1).*(-3*Prob1.M.c(i+1,1)+4*Prob1.M.c(i+1,2)-Prob1.M.c(i+1,3))/2/dxp.*...
        (-3*Prob1.M.Phi2(i+1,1)+4*Prob1.M.Phi2(i+1,2)-Prob1.M.Phi2(i+1,3))/2/dxp;
    
    Qohm_p(Np+2)= sig_eff(1)*((3*Prob1.M.Phi1(i+1,Np+2)-4*Prob1.M.Phi1(i+1,Np+1)+Prob1.M.Phi1(i+1,Np))/2/dxp).^2+Prob1.V.Keff_p1.*...
        ((3*Prob1.M.Phi2(i+1,Np+2)-4*Prob1.M.Phi2(i+1,Np+1)+Prob1.M.Phi2(i+1,Np))/2/dxp).^2 + ...
        2*Prob1.V.Keff_p1*R.*Prob1.M.T(i,Nal+Np+3)/F.*(1-Prob1.K.tplus(Np+2))./Prob1.M.c(i+1,Np+2).*(3*Prob1.M.c(i+1,Np+2)-4*Prob1.M.c(i+1,Np+1)+Prob1.M.c(i+1,Np))/2/dxp.*...
        (3*Prob1.M.Phi2(i+1,Np+2)-4*Prob1.M.Phi2(i+1,Np+1)+Prob1.M.Phi2(i+1,Np))/2/dxp;
    
    Qohm_n(1)   =  sig_eff(3)*((-3*Prob1.M.Phi1(i+1,Np+Ns+3)+4*Prob1.M.Phi1(i+1,Np+Ns+4)-Prob1.M.Phi1(i+1,Np+Ns+5))/2/dxn).^2+Prob1.V.Keff_n0.*...
        ((-3*Prob1.M.Phi2(i+1,Np+Ns+3)+4*Prob1.M.Phi2(i+1,Np+Ns+4)-Prob1.M.Phi2(i+1,Np+Ns+5))/2/dxn).^2 + ...
        2*Prob1.V.Keff_n0*R.*Prob1.M.T(i,Nal+Np+Ns+4)/F.*(1-Prob1.K.tplus(Np+Ns+3))./Prob1.M.c(i+1,Np+Ns+3).*(-3*Prob1.M.c(i+1,Np+Ns+3)+4*Prob1.M.c(i+1,Np+Ns+4)-Prob1.M.c(i+1,Np+Ns+5))/2/dxn.*...
        (-3*Prob1.M.Phi2(i+1,Np+Ns+3)+4*Prob1.M.Phi2(i+1,Np+Ns+4)-Prob1.M.Phi2(i+1,Np+Ns+5))/2/dxn;
    
    Qohm_n(Nn+2)= sig_eff(3)*((3*Prob1.M.Phi1(i+1,Np+Ns+Nn+4)-4*Prob1.M.Phi1(i+1,Np+Ns+Nn+3)+Prob1.M.Phi1(i+1,Np+Ns+Nn+2))/2/dxn).^2+Prob1.V.Keff_n1.*...
        ((3*Prob1.M.Phi2(i+1,Np+Ns+Nn+4)-4*Prob1.M.Phi2(i+1,Np+Ns+Nn+3)+Prob1.M.Phi2(i+1,Np+Ns+Nn+2))/2/dxn).^2 + ...
        2*Prob1.V.Keff_n1*R.*Prob1.M.T(i,Nal+Np+Ns+Nn+5)/F.*(1-Prob1.K.tplus(Np+Ns+Nn+4))./Prob1.M.c(i+1,Np+Ns+Nn+4).*(3*Prob1.M.c(i+1,Np+Ns+Nn+4)-4*Prob1.M.c(i+1,Np+Ns+Nn+3)+Prob1.M.c(i+1,Np+Ns+Nn+2))/2/dxn.*...
        (3*Prob1.M.Phi2(i+1,Np+Ns+Nn+4)-4*Prob1.M.Phi2(i+1,Np+Ns+Nn+3)+Prob1.M.Phi2(i+1,Np+Ns+Nn+2))/2/dxn;
    
    Qohm_s = Prob1.V.Keff_s.*((Prob1.M.Phi2(i+1,Is+1)-Prob1.M.Phi2(i+1,Is-1))/2/dxs).^2 + ...
        2*Prob1.V.Keff_s*R.*Prob1.M.T(i,Nal+Np+4:Nal+Np+Ns+3)/F.*(1-Prob1.K.tplus(Is))./Prob1.M.c(i+1,Is).*(Prob1.M.c(i+1,Is+1)-Prob1.M.c(i+1,Is-1))/2/dxs.*(Prob1.M.Phi2(i+1,Is+1)-Prob1.M.Phi2(i+1,Is-1))/2/dxs;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fTal(1)     = 0;
    fTal(2:Nal+1) = lambda_al/2/dxal^2*(Prob1.M.T(i,3:Nal+2)+Prob1.M.T(i,1:Nal)) + (rho_al*Cpal/dt-lambda_al/dxal^2)*Prob1.M.T(i,2:Nal+1) + Prob1.M.I(i+1)^2/sig_al;
    
    fTp(1) = 0;
    fTp(2:Np+1) = lambda_p/2/dxp^2*(Prob1.M.T(i,Nal+Ip+2)+Prob1.M.T(i,Nal+Ip)) + (rho_p*Cpp/dt-lambda_p/dxp^2)*Prob1.M.T(i,Nal+Ip+1) + Qrxn_p(2:Np+1) - Qrev_p(2:Np+1) + Qohm_p(2:Np+1);
    fTp(Np+2) = 0;
    
    fTs = lambda_s/2/dxs^2*(Prob1.M.T(i,Nal+Is+2)+Prob1.M.T(i,Nal+Is)) + (rho_s*Cps/dt-lambda_s/dxs^2)*Prob1.M.T(i,Nal+Is+1) + Qohm_s;
    
    fTn(1) = 0;
    fTn(2:Nn+1) = lambda_n/2/dxn^2*(Prob1.M.T(i,Nal+In+2)+Prob1.M.T(i,Nal+In)) + (rho_n*Cpn/dt-lambda_n/dxn^2)*Prob1.M.T(i,Nal+In+1) + Qrxn_n(2:Nn+1) - Qrev_n(2:Nn+1) + Qohm_n(2:Nn+1);
    fTn(Nn+2) = 0;
    
    fTco(1:Nco) = lambda_co/2/dxco^2*(Prob1.M.T(i,Nal+Np+Ns+Nn+7:Nal+Np+Ns+Nn+Nco+6)+Prob1.M.T(i,Nal+Np+Ns+Nn+5:Nal+Np+Ns+Nn+Nco+4)) + ...
        (rho_co*Cpco/dt-lambda_co/dxco^2)*Prob1.M.T(i,Nal+Np+Ns+Nn+6:Nal+Np+Ns+Nn+Nco+5) + Prob1.M.I(i+1)^2/sig_co;
    fTco(Nco+1) = 0;
    
    fT = [fTal fTp fTs fTn fTco]';
    Prob1.M.T(i+1,:) = iAT*fT;
    
end

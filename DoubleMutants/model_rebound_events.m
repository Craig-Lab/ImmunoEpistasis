function [time,sol,solstruc] = model_rebound_events(p,tspan,effects,prod_red,mut_struc)

%simulate until stabilized
opts = ddeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e-2);
sol1 = ddesd_f5(@ddefun,@(t,y) delayP(t,y,p),@history,[tspan(1) 200],opts);

%reset ICs
p.S0 = sol1.y(5,end);
p.S0_new = p.S0;
p.R0 = 0;
p.D0 = 0;
p.MPhi_R_0 = sol1.y(12,end);
p.MPhi_I_0 = sol1.y(13,end);
p.M0 = sol1.y(14,end);
p.N0 = sol1.y(15,end);
p.T0 = sol1.y(16,end);
p.L_U_0 = sol1.y(17,end);
p.L_B_0 = sol1.y(18,end);
p.G_U_0 = sol1.y(19,end);
p.G_B_0 = sol1.y(20,end);
p.C_U_0 = sol1.y(21,end);
p.C_B_0 = sol1.y(22,end);
p.F_U_0 = sol1.y(23,end);
p.F_B_0 = sol1.y(24,end);
p.A_0 = sol1.y(25,end);

%recalculate del_N to ensure S* = S(0)
p.del_N = (p.lam_S*(1-(p.S0+p.D0)/p.Smax)*p.S0)/(p.rho*p.S0*(p.N0^p.h_N/(p.N0^p.h_N+p.IC_50_N^p.h_N)));

solstruc = ddesd_f5(@ddefun,@(t,y) delayP(t,y,p),@history_virus,tspan,opts);
time = solstruc.x; %linspace(tspan(1),tspan(end),1000);
sol = solstruc.y; %deval(solstruc,time);

%------------------------------------------------------------------------
function [dydt] = ddefun(t,y,Z)

ylag1 = Z(:,1);     % tau_I eclipse phase lag
ylag2 = Z(:,2);     % tau_T lag for T-cell recruitment
ylag3 = Z(:,3);     % tau_R lag for resistant cell conversion back to susceptible cells
ylag4 = Z(:,3);     % tau_R lag for resistant cell conversion back to susceptible cells
       
V1 = y(1);%WT
V2 = y(2);%spike evasion strain
V3 = y(3);%IFN evasion strain
V4 = y(4);%double mutant
S = y(5);%susceptible cells
I1 = y(6);%WT infected
I2 = y(7);%spike evasion infected
I3 = y(8);%IFN evasion infected
I4 = y(9);%double mutant infected
R = y(10);%refractory cells
D = y(11);%damaged and dead cells
MPhi_R = y(12);%tissue-resident macrophages
MPhi_I = y(13);%inflammatory macrophages
M = y(14);%monocytes
N = y(15);%neutrophils
T = y(16);%CD8+ T cells
L_U = y(17);%unbound IL-6
L_B = y(18);%bound IL-6
G_U = y(19);%unbound GM-CSF
G_B = y(20);%bound GM-CSF
C_U = y(21);%unbound G-CSF
C_B = y(22);%bound G-CSF
F_U = y(23);%unbound Type I IFN
F_B = y(24);%bound Type I IFN
A = y(25);%neutralizing antibodies
C_BF=C_B/(p.A_C*N);%bound G-CSF fraction

TotalVL = V1+V2+V3+V4;
TotalInf = I1+I2+I3+I4;

%effects on production
lambda_F_prod = prod_red(1);
lambda_N_prod = prod_red(2);
lambda_M_prod = prod_red(3);
lambda_T_prod = prod_red(4);
lambda_A_prod = prod_red(5);
lambda_mono_prod = prod_red(6);
lambda_IL6_prod = prod_red(7);
lambda_GCSF_prod = prod_red(8);
lambda_GMCSF_prod = prod_red(9);

%mutation effects
gam_beta = effects(1);
gam_FI = effects(2);
gam_VN = effects(3);
gam_VM = effects(4);
gam_IN = effects(5);
gam_IM = effects(6);
gam_IT = effects(7);
gam_VA = effects(8);

Tkillrate = effects(:,7).*(p.del_I_T)/(p.K_del_T+TotalInf^p.h_T);

if TotalVL > 1e-5
    [dV1,~] = VirusRHS(p,I1,V1,V3,V2,V4,MPhi_I,N,A,effects(1,:));
    [dV2,~] = VirusRHS(p,I2,V2,V4,V1,V3,MPhi_I,N,A,effects(2,:));
    [dV3,~] = VirusRHS(p,I3,V3,V1,V4,V2,MPhi_I,N,A,effects(3,:));
    [dV4,~] = VirusRHS(p,I4,V4,V2,V3,V1,MPhi_I,N,A,effects(4,:));
else 
    V1 = 0; V2 = 0; V3 = 0; V4 = 0;
    dV1 = 0; dV2 = 0; dV3 = 0; dV4 = 0;
end

[dI1,b_I_w,m_I_w,inf_rate_w] = InfectedCellsRHS(p,I1,ylag1(5),ylag1(1),MPhi_I,N,T,F_B,Tkillrate(1),effects(1,:));
[dI2,b_I_m1,m_I_m1,inf_rate_m1] = InfectedCellsRHS(p,I2,ylag1(5),ylag1(2),MPhi_I,N,T,F_B,Tkillrate(2),effects(2,:));
[dI3,b_I_m2,m_I_m2,inf_rate_m2] = InfectedCellsRHS(p,I3,ylag1(5),ylag1(3),MPhi_I,N,T,F_B,Tkillrate(3),effects(3,:));
[dI4,b_I_m12,m_I_m12,inf_rate_m12] = InfectedCellsRHS(p,I4,ylag1(5),ylag1(4),MPhi_I,N,T,F_B,Tkillrate(4),effects(4,:));

dS = p.lam_S*(1-(S+TotalInf+D+R)/p.Smax)*S-(inf_rate_w*y(1)+inf_rate_m1*y(2)+inf_rate_m2*y(3)+inf_rate_m12*y(4)+p.rho*p.del_N*(N^p.h_N/(N^p.h_N+p.IC_50_N^p.h_N)))*S+p.d_R*ylag3(10);

%calculate reduced IFN effect
IFN_eff = (b_I_w*ylag1(1) + b_I_m1*ylag1(2) + b_I_m2*ylag1(3) + b_I_m12*ylag1(4));

dR = p.lam_S*(1-(S+TotalInf+D+R)/p.Smax)*R+(IFN_eff)-p.rho*p.del_N*(N^p.h_N/(N^p.h_N+p.IC_50_N^p.h_N))*R-p.d_R*ylag3(10);
dD = (m_I_w*I1+m_I_m1*I2+m_I_m2*I3+m_I_m12*I4)+(p.rho*p.del_N*(N^p.h_N/(N^p.h_N+p.IC_50_N^p.h_N)))*(S+R)-(p.d_D+(p.del_MPhi_D-p.del_D_MPhi)*(MPhi_R+MPhi_I))*D;

dMPhi_R = -p.a_I_MPhi*MPhi_R*(TotalInf+D)+(1-MPhi_R/p.MPhi_max)*p.lam_MPhi*MPhi_I/(TotalVL+p.eps_V_MPhi-p.del_MPhi_D*D*MPhi_R)-p.d_MPhi_R*MPhi_R;
dMPhi_I = lambda_M_prod*((p.p_MPhi_I_G*G_B^p.h_M_MPhi/(G_B^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)*M+p.p_MPhi_I_L*L_B/(L_B+p.eps_L_MPhi)*M))+p.a_I_MPhi*MPhi_R*(TotalInf+D)-p.del_MPhi_D*MPhi_I*D-(1-MPhi_R/p.MPhi_max)*p.lam_MPhi*MPhi_I/(TotalVL+p.eps_V_MPhi)-p.d_MPhi_I*MPhi_I;    
dM = lambda_mono_prod*((p.M_prod_star+(p.psi_M_max-p.M_prod_star)*G_B^p.h_M/(G_B^p.h_M+p.eps_G_M^p.h_M))*p.MR+p.p_M_I*TotalInf*M/(TotalInf+p.eps_I_M))-p.p_MPhi_I_G*G_B^p.h_M_MPhi*M/(G_B^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)-p.p_MPhi_I_L*L_B*M/(L_B+p.eps_L_MPhi)-p.d_M*M;

dN = lambda_N_prod*((p.N_prod_star+(p.psi_N_max-p.N_prod_star)*(C_BF-p.C_BF_star)/(C_BF-p.C_BF_star+p.eps_C_N))*p.NR+p.p_N_L*L_B/(L_B+p.eps_L_N))-p.d_N*N;

%calculate reduced T cell production
Tprod_new = p.p_T_I;
red_prod = Tprod_new'*sum(ylag2(6:9));

dT = lambda_T_prod*(red_prod/(1+L_B/p.eps_L_T)+(p.p_T_F*F_B/(F_B+p.eps_F_T))*T)-p.d_T*T;

dL_U = lambda_IL6_prod*(p.p_L_I*TotalInf/(TotalInf+p.eta_L_I)+p.p_L_MPhi*MPhi_I/(MPhi_I+p.eta_L_MPhi)+p.p_L_M*M/(M+p.eta_L_M))-p.k_lin_L*L_U-p.k_B_L*((N+T+M)*p.A_L-L_B)*L_U+p.k_U_L*L_B;
dL_B = -p.k_int_L*L_B+p.k_B_L*((N+T+ M)*p.A_L-L_B)*L_U-p.k_U_L*L_B;
dG_U = lambda_GMCSF_prod*(p.p_G_MPhi_I*MPhi_I/(MPhi_I+p.eta_G_MPhi)+p.p_G_M*M/(M+p.eta_G_M))-p.k_lin_G*G_U-p.k_B_G*(M*p.A_G-G_B)*G_U+p.k_U_G*G_B;
dG_B = -p.k_int_G*G_B+p.k_B_G*(M*p.A_G-G_B)*G_U-p.k_U_G*G_B;
dC_U = lambda_GCSF_prod*(p.p_C_M*M/(M+p.eta_C_M))-p.k_lin_C*C_U-p.k_B_C*(N*p.A_C-C_B)*C_U^p.stoch_C+p.k_U_C*C_B;
dC_B = -p.k_int_C*C_B+p.k_B_C*(N*p.A_C-C_B)*C_U^p.stoch_C-p.k_U_C*C_B;
dF_U = lambda_F_prod*(p.p_F_I*TotalInf/(TotalInf+p.eta_F_I)+p.p_F_MPhi*MPhi_I/(MPhi_I+p.eta_F_MPhi)+p.p_F_M*M/(M+p.eta_F_M))-p.k_lin_F*F_U-p.k_B_F*((T+S)*p.A_F-F_B)*F_U+p.k_U_F*F_B;
dF_B = -p.k_int_F*F_B+p.k_B_F*((T+S)*p.A_F-F_B)*F_U-p.k_U_F*F_B;

%calculate antibody evasion
p_A_new = lambda_A_prod*p.p_A*sum(ylag4(1:4));
neut = ((effects(:,8)*p.del_V_A)*(A).^p.h_A)./(p.eps_V_A^p.h_A+(A).^p.h_A);
neut_new = neut'*(y(1:4));

dA = p_A_new-p.d_A*A-neut_new;

dydt = [dV1;dV2;dV3;dV4;dS;dI1;dI2;dI3;dI4;dR;dD;dMPhi_R;dMPhi_I;dM;dN;dT;dL_U;dL_B;dG_U;dG_B;dC_U;dC_B;dF_U;dF_B;dA];

end

%------------------------------------------------------------------------
function s = history(t)
         s = [0;0;0;0;p.S0;0;0;0;0;p.R0;p.D0;p.MPhi_R_0;p.MPhi_I_0;p.M0;p.N0;p.T0;p.L_U_0;p.L_B_0;p.G_U_0;p.G_B_0;p.C_U_0;p.C_B_0;p.F_U_0;p.F_B_0;p.A_0];
end

function s = history_virus(t)
    if t < 0
          s = [0;0;0;0;p.S0;0;0;0;0;p.R0;p.D0;p.MPhi_R_0;p.MPhi_I_0;p.M0;p.N0;p.T0;p.L_U_0;p.L_B_0;p.G_U_0;p.G_B_0;p.C_U_0;p.C_B_0;p.F_U_0;p.F_B_0;p.A_0];

    else % Initial values at t >= 0
          s = [p.V10;p.V20;p.V30;p.V40;p.S0;p.I10;p.I20;p.I30;p.I40;p.R0;p.D0;p.MPhi_R_0;p.MPhi_I_0;p.M0;p.N0;p.T0;p.L_U_0;p.L_B_0;p.G_U_0;p.G_B_0;p.C_U_0;p.C_B_0;p.F_U_0;p.F_B_0;p.A_0];
    end
end

function d = delayP(t,y,p)
%This function sets up the delay vectors necessary for the DDE solver.
d = [t-p.tau_I,t-p.tau_T,t-p.tau_R,t-p.tau_A];     
end
%-------------------------------------------------------------------------
end
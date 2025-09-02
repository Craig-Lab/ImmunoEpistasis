function [Healthy,Immunodeficiency] = RunEpistasisModel_double_mut(V10,V20,V30,V40,I10,I20,I30,I40,cGamma,mu_F,mu_I,mu_D,effects,prod_effects,mut_struc)

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% Healthy individuals

p = load_parameters();

%virus specific parameters
%initial conditions
p.V10 = V10;%WT
p.V20 = V20;%spike evasion strain
p.V30 = V30;%immune evasion strain
p.V40 = V40;%double mutation

p.I10 = I10;%WT
p.I20 = I20;%spike evasion strain
p.I30 = I30;%immune evasion strain
p.I40 = I40;%double mutation

%multistrain evolutionary model parameters
p.cGamma = cGamma;%mutation cost
p.mu_F = mu_F;%rate of mutation to immune evasive strain
p.mu_I = mu_I;%rate of mutation to spike evasive strain
p.mu_D = mu_D;%rate of double mutations

%---------------------------------------------------

%refractory cell entry (if not uncommented, d_R > 0)
%p.d_R = 0;

%set immune parameters
prod_red = ones(size(prod_effects,2),1);

%-----------------------------------------------------------------------

%calculate homeostasis values and check for issues
p = Homeostasis_calculations(p);

estimated_params = [p.p_L_M  p.L_B_star  p.MPhi_I_star  p.p_L_MPhi  p.eta_G_MPhi  p.p_G_MPhi_I  p.T_prod_star p.T_M_prod_star p.p_G_M...
    p.M_prod_star p.N_prod_star p.eta_F_MPhi];

if isempty(find(estimated_params<0))==0
    disp('Negative parameter')
elseif isempty(find(estimated_params>1e9))==0
    disp('Extremely large parameter')
end

%----------------------------------------------------------------------

%solve model
tspan = [0 150];
disp('Healthy')

IC1 =[p.V10;p.V20;p.V30;p.V40;p.S0;p.I10;p.I20;p.I30;p.I40;p.R0;p.D0;p.MPhi_R_0;p.MPhi_I_0;p.M0;p.N0;p.T0;p.L_U_0;p.L_B_0;p.G_U_0;p.G_B_0;...
    p.C_U_0;p.C_B_0;p.F_U_0;p.F_B_0;p.A_0;0];

[time1,sol1,solstruc1] = model_rebound_events(p,tspan,effects,prod_red,mut_struc);
time_deval = linspace(tspan(1),tspan(2),1e3);
sol1_deval = interp1(time1,sol1',time_deval);  

%save results
Healthy.time_deval = time_deval;
Healthy.sol_deval = sol1_deval;
Healthy.IC = IC1;
Healthy.p = p;

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% Immunodeficiency
p = load_parameters();

%---------------------------------------------------

%virus specific parameters
%initial conditions
p.V10 = V10;%WT
p.V20 = V20;%spike evasion strain
p.V30 = V30;%immune evasion strain
p.V40 = V40;%double mutation

p.I10 = I10;%WT
p.I20 = I20;%spike evasion strain
p.I30 = I30;%immune evasion strain
p.I40 = I40;%double mutation

%multistrain evolutionary model parameters
p.cGamma = cGamma;%mutation cost
p.mu_F = mu_F;%rate of mutation to immune evasive strain
p.mu_I = mu_I;%rate of mutation to spike evasive strain
p.mu_D = mu_D;%rate of double mutations

%-----------------------------------------------------------------------

%refractory cell entry (if not uncommented, d_R > 0)
%p.d_R = 0;

%set immune parameters
prod_red = prod_effects;

%-----------------------------------------------------------------------

%calculate homeostasis values and check for issues

%p = Homeostasis_calculations(p);

lambda_F_prod = prod_red(1);
lambda_N_prod = prod_red(2);
lambda_M_prod = prod_red(3);
lambda_T_prod = prod_red(4);
lambda_A_prod = prod_red(5);
lambda_mono_prod = prod_red(6);
lambda_IL6_prod = prod_red(7);
lambda_GCSF_prod = prod_red(8);
lambda_GMCSF_prod = prod_red(9);

p.MPhi_R_0 = lambda_M_prod*p.MPhi_R_0;
p.MPhi_I_0 = lambda_M_prod*p.MPhi_I_0;
p.N0 = lambda_N_prod*p.N0;
p.T0 = lambda_T_prod*p.T0;
p.F_U_0 = lambda_F_prod*p.F_U_0;
p.F_B_0 = lambda_F_prod*p.F_B_0;
p.A_0 = lambda_A_prod*p.A_0;
p.M0 = lambda_mono_prod*p.M0;
p.L_U_0 = lambda_IL6_prod*p.L_U_0;
p.L_B_0 = lambda_IL6_prod*p.L_B_0;
p.C_U_0 = lambda_GCSF_prod*p.C_U_0;
p.C_B_0 = lambda_GCSF_prod*p.C_B_0;
p.G_U_0 = lambda_GMCSF_prod*p.G_U_0;
p.G_B_0 = lambda_GMCSF_prod*p.G_B_0;

p = Homeostasis_calculations(p);

if lambda_M_prod ~= 1
    p.MPhi_I_0 = lambda_M_prod*p.MPhi_I_0;
end

if lambda_GMCSF_prod ~= 1
    p.G_U_0 = lambda_GMCSF_prod*p.G_U_0;
    p.G_B_0 = lambda_GMCSF_prod*p.G_B_0;
    p.G_U_star = p.G_U_0;    % homeostatic unbound concentration
    p.G_B_star = (p.k_B_G*p.M_star*p.A_G*p.G_U_star)/(p.k_int_G+p.k_U_G+p.k_B_G*p.G_U_star);  % homeostatic bound concentration
    p.G_B_0 = p.G_B_star;       % Initial G-CSF (bound) in ng/mL
    p.MPhi_I_star = (p.p_MPhi_I_G*p.G_B_star^p.h_M_MPhi*p.M_star/(p.G_B_star^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)...
                    +p.p_MPhi_I_L*p.L_B_star*p.M_star/(p.L_B_star+p.eps_L_MPhi))/((1-p.MPhi_R_star/p.MPhi_max)*p.lam_MPhi/p.eps_V_MPhi+p.d_MPhi_I);
    p.MPhi_I_0 = p.MPhi_I_star;
    p.eta_G_MPhi = p.eta_L_MPhi;
    p.p_G_MPhi_I = -(-p.k_lin_G*p.G_U_star-p.k_B_G*(p.M_star*p.A_G-p.G_B_star)*p.G_U_star+p.k_U_G*p.G_B_star)/(p.MPhi_I_star/(p.MPhi_I_star+p.eta_G_MPhi)...
                   +p.M_star/(p.M_star+p.eta_G_M));
    p.p_G_M = p.p_G_MPhi_I; %this was assuming that p_G_MPhi_I = p_G_M
    p.p_C_M = p.p_G_M/100; %unknown 
    DD = 1/p.MR*(p.p_MPhi_I_G*p.G_B_star^p.h_M_MPhi*p.M_star/(p.G_B_star^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)...
        +p.p_MPhi_I_L*p.L_B_star*p.M_star/(p.L_B_star+p.eps_L_MPhi)+p.d_M*p.M_star);
    EE = p.G_B_star^p.h_M/(p.G_B_star^p.h_M+p.eps_G_M^p.h_M);
    p.M_prod_star = (DD-p.psi_M_max*EE)/(1-EE);
end

IC2 =[p.V10;p.V20;p.V30;p.V40;p.S0;p.I10;p.I20;p.I30;p.I40;p.R0;p.D0;p.MPhi_R_0;p.MPhi_I_0;p.M0;p.N0;p.T0;p.L_U_0;p.L_B_0;p.G_U_0;p.G_B_0;...
    p.C_U_0;p.C_B_0;p.F_U_0;p.F_B_0;p.A_0;0];

estimated_params = [p.p_L_M  p.L_B_star  p.MPhi_I_star  p.p_L_MPhi  p.eta_G_MPhi  p.p_G_MPhi_I  p.T_prod_star p.T_M_prod_star p.p_G_M...
    p.M_prod_star p.N_prod_star p.eta_F_MPhi];

if isempty(find(estimated_params<0))==0
    disp('Negative parameter')
elseif isempty(find(estimated_params>1e9))==0
    disp('Extremely large parameter')
end

%----------------------------------------------------------------------

%solve model
%tspan = [0 300];
tspan = [0 150];
disp('Immunodeficiency')

[time2,sol2,solstruc2] = model_rebound_events(p,tspan,effects,prod_red,mut_struc);

time_deval = linspace(tspan(1),tspan(2),1e3);
sol2_deval = interp1(time2,sol2',time_deval);

%save results
Immunodeficiency.time_deval = time_deval;
Immunodeficiency.sol_deval = sol2_deval;
Immunodeficiency.IC = IC2;
Immunodeficiency.p = p;

function [RHS_Virus,m_V] = VirusRHS(p,I,V,Vmut1,Vmut2,Vmut3,MPhi_I,N,A,effects)

%mutation effects
gam_beta = effects(1);
gam_FI = effects(2);
gam_VN = effects(3);
gam_VM = effects(4);
gam_IN = effects(5);
gam_IM = effects(6);
gam_IT = effects(7);
gam_VA = effects(8);

m_V = (gam_VM*p.del_V_MPhi).*(MPhi_I) + (gam_VN*p.del_V_N).*(N) + p.d_V_spec...
    + ((gam_VA*p.del_V_A).*(A).^p.h_A)./(p.eps_V_A.^p.h_A + (A).^p.h_A);

RHS_Virus = p.phat*I - m_V.*V  + p.mu_F*(Vmut1 - V) + p.mu_I*(Vmut2 - V) + p.mu_D*(Vmut3 - V);% 

end
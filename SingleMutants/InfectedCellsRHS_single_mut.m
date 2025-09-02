function [RHS_InfectedCells,b_I,m_I,inf_rate] = InfectedCellsRHS_single_mut(p,I,Slag,Vlag,MPhi_I,N,T,F_B,Tkillrate,lambda,effects)

%[RHS_InfectedCells,b_I,m_I,inf_rate] = InfectedCellsRHS_single_mut(p,I,Imut1,Imut2,Imut3,Slag,Vlag,MPhi_I,N,T,F_B,Tkillrate,lambda,effects)

%mutation effects
gam_beta = effects(1);
gam_FI = effects(2);
gam_VN = effects(3);
gam_VM = effects(4);
gam_IN = effects(5);
gam_IM = effects(6);
gam_IT = effects(7);
gam_VA = effects(8);

%immunodeficiency values
lambda_F = lambda(1);
lambda_N = lambda(2);
lambda_M = lambda(3);
lambda_T = lambda(4);
lambda_A = lambda(5);

inf_rate = (1/gam_beta)*p.beta;
b_I = ((1/gam_beta)*p.beta).*((((1/gam_FI)*p.eps_F_I))./(((1/gam_FI)*p.eps_F_I)+(F_B))).*Slag;
m_I = p.d_I+gam_IN*((p.del_N)./(1+(p.IC_50_N./(N)).^p.h_N))+gam_IM*(p.del_I_MPhi.*MPhi_I)+Tkillrate.*T;

RHS_InfectedCells = b_I.*Vlag - m_I.*I;

end
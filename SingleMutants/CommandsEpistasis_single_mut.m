% Main commands code to run within-host evolution model (single mutant)

clearvars
close all

production_effect = 1;%if 1 -> effect is active, if 0 -> effect is inactive

% changes for production deficiencies
% lambda1 -> IFN effect 
% lambda2 -> neutrophil effect 
% lambda3 -> macrophage effect 
% lambda4 -> T cell effect 
% lambda5 -> antibody effect
% lambda6 -> monocyte effect
% lambda7 -> IL-6 effect
% lambda8 -> G-CSF effect
% lambda9 -> GM-CSF effect

for jj = 1:9

    def_red = 1;

    if production_effect == 1 && jj == 1
        prod_vals = 0.65*def_red;
    elseif production_effect == 1 && jj == 3
        prod_vals = 0.65*def_red;
    elseif production_effect == 1 && jj == 4
        prod_vals = 0.25*def_red;
    elseif production_effect == 1 && jj == 5
        prod_vals = 0.1*def_red;
    elseif production_effect == 1 && jj == 6
        prod_vals = 0.65*def_red;
    elseif production_effect == 1 && jj == 9
        prod_vals = 0.65*def_red;
    elseif production_effect == 1
        prod_vals = 0.65*def_red;
    elseif production_effect == 0
        prod_vals = 1;
    end

    prod_effects = ones(1,9);
    prod_effects(jj) = prod_vals;

    lambda_vals = [1,1,1,1,1];%keep as 1s
    immune_names = {'IFN_deficiency','Neutrophil_deficiency','Macrophage_deficiency','T_cell_deficiency','Antibody_deficiency'};
    
    for nn = 1
    
        %initialize immunodeficiencies vector
        immune_effects = ones(size(lambda_vals,2),1); 
        %set immunodeficiency
        immune_effects(nn) = lambda_vals(nn);
    
        %---------------------------------------------------
       
        % effects on viral evasion of immune system (all values between 0 and 1)
        % gamma1 -> beta
        % gamma2 -> epsilonFI
        % gamma3 -> deltaVN
        % gamma4 -> deltaVM
        % gamma5 -> deltaIN
        % gamma6 -> deltaIM
        % gamma7 -> deltaIT
        % gamma8 -> deltaVA

        strength_selection_vals = ones(4,8);
        gam_val = [0.95,0.5,0.6,0,0,0,0.001,0.4];

        for mm = [1,2,3,7,8]%1:size(strength_selection_vals,2) 
            
            disp(jj)
            disp(nn)
            disp(mm)
    
            strength_selection = strength_selection_vals(mm);
            effects = strength_selection_vals;
            effects(2,mm) = gam_val(mm);
    
            %total viral load and proportions for WT and mutant
            V0 = 4.5;
            init_freq = 1e-2;

            %initial conditions
            Mut0 = V0*init_freq;
            WT0 = V0*(1-init_freq);
            DoubleMut0 = V0*((init_freq^2)/(1-init_freq));
            
            V10 = WT0;%wild type 
            V20 = Mut0;%mutant 1
            V30 = 0;%mutant 2
            V40 = 0;%double mutant
            
            %infected cells
            I10 = 0;%WT strain
            I20 = 0;%strain with higher infectivity
            I30 = 0;%strain with immune escape
            I40 = 0;%double mutant
    
            %multistrain evolutionary model parameters
            cGamma = 0;%cost of mutation
            mu_F = 0;%rate of mutation to immune evasion strain
            mu_I = 0;%rate of mutation to strain with higher infectivity
            mu_D = 0;%rate of double mutations
            
            %---------------------------------------------------
            %run model
            [Healthy,Immunodeficiency] = RunEpistasisModel_single_mut(immune_effects,V10,V20,V30,V40,I10,I20,I30,I40,cGamma,mu_F,mu_I,mu_D,effects,prod_effects);
       
        end
    
    end

end
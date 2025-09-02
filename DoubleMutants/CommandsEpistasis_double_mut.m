% Main commands code to run within-host evolution model (double mutants)

clearvars
close all

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

production_effect = 1;%if 1 -> effect is active, if 0 -> effect is inactive
for jj = 1:9

    if production_effect == 1 && jj == 1
        prod_vals = 0.65;
    elseif production_effect == 1 && jj == 3
        prod_vals = 0.65;
    elseif production_effect == 1 && jj == 4
        prod_vals = 0.25;
    elseif production_effect == 1 && jj == 5
        prod_vals = 0.1;
    elseif production_effect == 1 && jj == 6
        prod_vals = 0.65;
    elseif production_effect == 1 && jj == 9
        prod_vals = 0.65;
    elseif production_effect == 1
        prod_vals = 0.65;
    elseif production_effect == 0
        prod_vals = 1;
    end

    prod_effects = ones(1,9);
    prod_effects(jj) = prod_vals;

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
    mutation_names = {'\beta','\epsilon_{FI}','\delta_{VN}','\delta_{VM}','\delta_{IN}','\delta_{IM}','\delta_{IT}','\delta_{VA}'};
    gam_val = [0.7,0.5,0,0,0,0,0,0];

    %first mutation
    for mm = 1:size(strength_selection_vals-1,2)

        if gam_val(mm) == 0%skip weak mutations
            continue
        end
            
        %second mutation
        for nn = 1:size(strength_selection_vals,2)
    
            if gam_val(nn) == 0%skip weak mutations
                continue
            end

            disp([mutation_names{mm},mutation_names{nn}])
            strength_selection = strength_selection_vals(mm);

            %initial viral loads
            V0 = 4.5;
            init_freq = 1e-2;

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

            for muts = 1:3

                %mut_struc = 1 -> just first mutation
                %mut_struc = 2 -> just second mutation
                %mut_struc = 3 -> both mutations

                mut_struc = muts;
    
                if mut_struc == 1
    
                    effects = strength_selection_vals;
                    effects(2,mm) = gam_val(mm);%first mutant
                    effects(4,mm) = gam_val(mm);%double mutant 1
    
                    V10 = V0*((1-init_freq));%WT 
                    V20 = V0*(init_freq);%mutant 1
                    V30 = 0;%mutant 2
                    V40 = 0;%double mutant

                    %run model
                    [Healthy10,Immunodeficiency10] = RunEpistasisModel_double_mut(V10,V20,V30,V40,I10,I20,I30,I40,cGamma,mu_F,mu_I,mu_D,effects,prod_effects,mut_struc);
                
                elseif mut_struc == 2
    
                    effects = strength_selection_vals;
                    effects(3,nn) = gam_val(nn);%second mutant
                    effects(4,nn) = gam_val(nn);%double mutant 2
                    
                    V10 = V0*((1-init_freq));%WT
                    V20 = 0;%mutant 1
                    V30 = V0*(init_freq);%mutant 2
                    V40 = 0;%double mutant

                    %run model
                    [Healthy02,Immunodeficiency02] = RunEpistasisModel_double_mut(V10,V20,V30,V40,I10,I20,I30,I40,cGamma,mu_F,mu_I,mu_D,effects,prod_effects,mut_struc);

                elseif mut_struc == 3
    
                    effects = strength_selection_vals;
                    effects(2,mm) = gam_val(mm);%first mutant
                    effects(3,nn) = gam_val(nn);%second mutant
                    effects(4,mm) = gam_val(mm);%double mutant 1
                    effects(4,nn) = gam_val(nn);%double mutant 2
    
                    V10 = V0*((1-init_freq)^2);%WT
                    V20 = V0*(init_freq*(1-init_freq));%mutant 1
                    V30 = V0*(init_freq*(1-init_freq));%mutant 2
                    V40 = V0*(init_freq^2);%double mutant

                    %run model
                    [Healthy12,Immunodeficiency12] = RunEpistasisModel_double_mut(V10,V20,V30,V40,I10,I20,I30,I40,cGamma,mu_F,mu_I,mu_D,effects,prod_effects,mut_struc);
                end
            end 
        end
   end
end
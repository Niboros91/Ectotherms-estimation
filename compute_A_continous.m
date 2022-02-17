function A_tot = compute_A_continous(Pest_stages)
    %Input =  Parameters for each stage
    %Output =  Continuous matrix A
    
    % We have 9 states
    A_tot = [-Pest_stages(1).growth-Pest_stages(1).death 0 0 0 0 0 0 Pest_stages(8).birth;... %Egg (1)
        Pest_stages(1).growth -Pest_stages(2).growth-Pest_stages(2).death 0 0 0 0 0 0;... %L1 (2)
        0 Pest_stages(2).growth -Pest_stages(3).growth-Pest_stages(3).death 0 0 0 0 0;... %L2 (3)
        0 0 Pest_stages(3).growth -Pest_stages(4).growth-Pest_stages(4).death 0 0 0 0;... %L3 (4)
        0 0 0 Pest_stages(4).growth -Pest_stages(5).growth-Pest_stages(5).death 0 0 0;...%P (5)
        0 0 0 0 Pest_stages(5).sex_ratio*Pest_stages(6).growth -Pest_stages(6).growth-Pest_stages(6).death 0 0;...%AM (6)
        0 0 0 0 (1-Pest_stages(6).sex_ratio)*Pest_stages(7).growth 0 -(Pest_stages(7).mate-Pest_stages(7).death) 0;%NMF (7)
        0 0 0 0 0 0 Pest_stages(7).mate-Pest_stages(7).death -Pest_stages(8).growth-Pest_stages(8).death]; %MF (8)

    
end
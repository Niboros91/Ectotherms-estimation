function Pest_stages = Initialize_stages_ode(Birth,Death,Growth,S_R,R_mate,R_remate,Pest_stages)
    %Input =  Different rates at each instant
    %Output =  Parameters associated to each stage

    for i=1:length(Pest_stages) %All stages have the same parameters
        Pest_stages(i).birth = Birth;
        Pest_stages(i).death = Death;
        Pest_stages(i).growth = Growth;
        Pest_stages(i).sex_ratio = S_R;
        Pest_stages(i).mate = R_mate;
        Pest_stages(i).remate = R_remate;

    end

end
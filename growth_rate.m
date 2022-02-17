function R = growth_rate(T,a,T_l,T_m,m)
    %Input =  fitting parameters and temperature T
    %Output =  rate for the given temperature
    
    R = a*T*(T-T_l)*((T_m-T)^(1/m));

end
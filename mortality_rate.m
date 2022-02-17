function M = mortality_rate(T,a1,b1,c1,d1,e1)
    %Input =  fitting parameters and temperature T
    %Output =  rate for the given temperature
    
    % Data from Ryan et al. (2016) fitted with a fourth polynomial function
    % (script fitmort.py)
    
    if T <=6 && T >= 32 %The function is bounded
        M = 1;
    else
        M= a1 * T^4 + b1 * T^3 + c1 * T^2 + d1 * T + e1;

    end
end
    
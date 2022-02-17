function B = birth_rate_suzuki(T,alpha,gamma,lambda,delta,tau)
    %Input =  fitting parameters and temperature T
    %Output =  rate for the given temperature
    
    % Data and function from Ryan et al. (2016)
    Tmin = 5; %Minimum temperature for oviposition
    Tmax = 30; %Maximum temperature for oviposition
    
    if T <= Tmin && T >= Tmax
        B = 0;
    else
    B = alpha * ((gamma + 1)/(pi * lambda^(2 * gamma + 2)) * (lambda^2 - ((T - tau)^2 + delta^2))^gamma);
    end
end
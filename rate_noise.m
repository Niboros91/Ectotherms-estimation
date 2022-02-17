function [w_d, w_m, w_f] = rate_noise(T,v_a1,v_b1,v_c1,v_d1,v_e1,v_a,v_T_l,v_T_m,v_m,a,T_l,T_m,m)

%Mortality rate error
w_m= ((v_a1 * T^4)^2 + (v_b1 * T^3)^2 + (v_c1 * T^2)^2 + (v_d1 * T)^2 + (v_e1)^2)^(1/2);

%Growth rate error
R = a*T*(T-T_l)*((T_m-T)^(1/m));

w_d = (((R/a) * v_a)^2  + ((R/(T-T_l)) * v_T_l)^2 + ( (R/m)*(T_m - T)^(1/(m*(m-1)))  *v_T_m  )^2 + ( (R/m^2) * log(T_m -T)*v_m )^2)^(1/2);

%Fertility rate error (randomly generated based on the other two errors)
w_f = min(w_m,w_d) + (max(w_m,w_d)-min(w_m,w_d)) .* rand(1,1);
end
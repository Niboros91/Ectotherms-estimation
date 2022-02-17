% ODE model for Drosophila suzukii
%Initialization
clear, clc % To clear the workspace and the command function

%% PARAMETERS functions

%Growth rate (Drosophila suzukii)
 a = 1.2*(10^(-4));
 T_l = 3;
 T_m = 30;
 m = 6;
 
%Uncertainty growth rate
v_a = 0.15 * (10^(-4));
v_T_l = 2;
v_T_m = 1;
v_m = 3;

%Mortality rate (Drosophila suzukii)

a1 = -5.4E-06;
b1 = 0.0005194;
c1 = -0.0116827;
d1 = 2.16E-05;
e1 = 1.3146586;

%Uncertainty mortality
v_a1 = 1.4E-06;
v_b1 = 0.0008;
v_c1 = 0.02;
v_d1 = 0.3E-05;
v_e1 = 0.9;

% Birth rate (Drosophila suzukii)
alpha = 659.06;
gamma = 88.53;
lambda = 52.32;
delta = 6.06;
tau = 22.87;

% Sex ratio
S_R = 0.5;

% Mating ratio
R_remate = 0;
R_mate = 1;

% Additional parameters
N_stages = 8; % Eggs/3 larva stages/ pupa /male/ unmated female/ mated female

%% ITERATIONS

Number_iterations = 2;

for it = 1:Number_iterations
    %% DATA

    %Temperature (Obtained from PANTHEON 2020)
    load('data_pantheon_2020'); %Temperature is Temp_avg
    temperature(277:288) = 15; %Temperature from 277 to 288 missing (malfunctioning of weather station)
    Temp_avg = temperature(92:end); %Selecting data fromt the 1st of April
    
    
    %% NOISE

    % Temperature measurement
    %Small temperature error in the measurement
    max_t = 0.2;
    min_t = -0.2;
    t_var = min_t + (max_t-min_t) .* rand(length(Temp_avg),1);

    %We separate between the temperature measured and the real temperature
    Temp_meas = Temp_avg;
    Temp_avg = Temp_avg + t_var';


    %% RATES COMPUTATION


    % REAL RATES
    %Compute the rates based on temperature
    G_R = growth_rate(Temp_avg(1),a,T_l,T_m,m); %Calling the growth rate function

    B_R = birth_rate_suzuki(Temp_avg(1),alpha,gamma,lambda,delta,tau);%Calling the birth rate function

    M_R = mortality_rate(Temp_avg(1),a1,b1,c1,d1,e1); %Calling the mortality rate function

    %Initialize stages (synthetic data)
    Pest_stages(1:N_stages) = stage; % We create an array of the stage class
    Pest_stages = Initialize_stages_ode(B_R,M_R,G_R,S_R,R_mate,R_remate,Pest_stages); %We initialize the parameters associated to each stage class



    % MEASURED RATES
    %Compute the rates based on temperature
    G_R_op = growth_rate(Temp_meas(1),a,T_l,T_m,m); %Calling the growth rate function

    B_R_op = birth_rate_suzuki(Temp_meas(1),alpha,gamma,lambda,delta,tau);%Calling the birth rate function

    M_R_op = mortality_rate(Temp_meas(1),a1,b1,c1,d1,e1); %Calling the mortality rate function
    
    [w_d, w_m, w_f] = rate_noise(Temp_avg(1),v_a1,v_b1,v_c1,v_d1,v_e1,v_a,v_T_l,v_T_m,v_m,a,T_l,T_m,m);
    
    
    % Development noise
    max_meas = w_d;
    min_meas = -w_d;
    G_R_noise = min_meas + (max_meas-min_meas) .* rand(1,1);
    G_R_op = G_R_op + G_R_noise;
    
    % Birth noise
    max_meas = w_f;
    min_meas = -w_f;
    B_R_noise = min_meas + (max_meas-min_meas) .* rand(1,1);
    B_R_op = B_R_op + B_R_noise;
    
    % Mortality noise
    max_meas = w_m;
    min_meas = -w_m;
    M_R_noise = min_meas + (max_meas-min_meas) .* rand(1,1);
    M_R_op = M_R_op + M_R_noise;
    

    %Initialize stages (estimation)
    Pest_stages_op(1:N_stages) = stage; % We create an array of the stage class
    Pest_stages_op = Initialize_stages_ode(B_R_op,M_R_op,G_R_op,S_R,R_mate,R_remate,Pest_stages_op); %We initialize the parameters associated to each stage class



    %% DYNAMIC MODEL

    %Initial values for the variables
    x = zeros(N_stages,1); %Array of variables

    %Real matrix
    A_cont =compute_A_continous(Pest_stages); %We build the continuous time A matrix
    % Discretization of the system based on the zero order holder
    sysc = ss(A_cont,[],eye(8),[]);
    sysd = c2d(sysc,1,'zoh');
    A_dis = sysd.A; %Discretized matrix A

    %Open loop matrix
    A_cont_op =compute_A_continous(Pest_stages_op); %We build the continuous time A matrix
    % Discretization of the system based on the zero order holder
    sysc = ss(A_cont_op,[],eye(8),[]);
    sysd = c2d(sysc,1,'zoh');
    A_dis_op = sysd.A; %Discretized matrix A

    %% MEASUREMENT MODEL

    Trap_eff = 0.5;

    % We assume that we can measure the 3 adult states

    C = [zeros(1,5), Trap_eff, 0, 0; zeros(1,6), Trap_eff, 0; zeros(1,7), Trap_eff];
    
    % We have a noise in the measurement
    max_meas = 2; %upper bound noise
    min_meas = -2; %lower bound noise
    meas_noise = min_meas + (max_meas-min_meas) .* rand(3,length(Temp_avg));


    %% MEASUREMENT PROCESS

    Number_meas = 20; %We do 15 random measurements

    %Certain days where there are measurements
    test = randperm(length(Temp_avg));
    measures = sort(test(1:Number_meas));

    counter =1; %Initialize counter of measurements
    %% SIMULATIONS
    Simulation_time =length(Temp_avg); %Simulation lenght based on the temperature array introduced

    %Initial conditions
    x(1) = 1000000; % eggs
    x(7) = 0; % adult males
    x(8) = 1000000; % adult mated females
    
    %Variation in the initial conditions
    max_x = 10000;
    min_x = -10000;
    x_var = min_x + (max_x-min_x) .* rand(2,1);

    %Initial estimation (Open loop and both EKFs)
    x_open = x;
    x_open(1) = 1000000 + x_var(1); % eggs
    x_open(8) = 1000000 + x_var(2); % adult mated females
    
    x_est = x;
    x_est(1) = 1000000 + x_var(1); % eggs
    x_est(8) = 1000000 + x_var(2); % adult mated females
    
    %Initial conditions for the estimator
    x_ini =x_est;
    
    % Historic evolution of the states
    x_hist= x; %We start the array x_hist to keep a historic of the evolution of the states
    x_open_hist = x_open;
    x_est_hist = x_est;
    
    
   L = [-x_est(1) -x_est(1) x_est(8) zeros(1,12) ;... %Egg (1)
    x_est(1) 0 0 -x_est(2) -x_est(2) zeros(1,10);... %L1 (2)
    zeros(1,3) x_est(2) 0 -x_est(3) -x_est(3) zeros(1,8);... %L2 (3)
    zeros(1,5) x_est(3) 0 -x_est(4) -x_est(4) zeros(1,6);... %L3 (4)
    zeros(1,7) x_est(4) 0 -x_est(5) -x_est(5) zeros(1,4);...%P (5)
    zeros(1,9) 0.5*x_est(5) 0 -x_est(6) zeros(1,3);...%AM (6)
    zeros(1,9) 0.5*x_est(5) 0 0  -x_est(7) -x_est(7) 0;%NMF (7)
    zeros(1,12) x_est(7) 0 -x_est(8)]; %MF (8)
    

    
    for t=1:(Simulation_time-1) %Loop for the simulation

        x = A_dis*x; %Compute the state at the next time step
        x_open = A_dis_op*x_open; %Compute the state at the next time step

        measuring = 0;
        y=0;
        if counter <= length(measures) %Measurements only arrive at certain instants
            if measures(counter) == t
                counter = counter +1;
                y = C*x + meas_noise(t+1);
                measuring = 1;
            end
        end
        x_est = EKF_noise(t,A_dis_op,L,C,N_stages,x_ini,y,measuring);

        x_hist = [x_hist,x]; %Store the states
        x_open_hist = [x_open_hist,x_open]; %Store the states
        x_est_hist = [x_est_hist,x_est]; %Store the states

        % Recompute rates for the new conditions
        G_R = growth_rate(Temp_avg(t+1),a,T_l,T_m,m);
        B_R = birth_rate_suzuki(Temp_avg(t+1),alpha,gamma,lambda,delta,tau);
        M_R = mortality_rate(Temp_avg(t+1),a1,b1,c1,d1,e1);
        Pest_stages = Initialize_stages_ode(B_R,M_R,G_R,S_R,R_mate,R_remate,Pest_stages);
        
        % Estimate rates from measurements
        G_R_op = growth_rate(Temp_meas(t+1),a,T_l,T_m,m); %Calling the growth rate function
        B_R_op = birth_rate_suzuki(Temp_meas(t+1),alpha,gamma,lambda,delta,tau);%Calling the birth rate function
        M_R_op = mortality_rate(Temp_meas(t+1),a1,b1,c1,d1,e1); %Calling the mortality rate function
        Pest_stages_op = Initialize_stages_ode(B_R_op,M_R_op,G_R_op,S_R,R_mate,R_remate,Pest_stages_op); %We initialize the parameters associated to each stage class
        
        %Rates uncertainty based on current conditions
        [w_d, w_m, w_f] = rate_noise(Temp_avg(t+1),v_a1,v_b1,v_c1,v_d1,v_e1,v_a,v_T_l,v_T_m,v_m,a,T_l,T_m,m);
       
        % Development noise
        max_meas = w_d;
        min_meas = -w_d;
        G_R_noise = min_meas + (max_meas-min_meas) .* rand(1,1);
        G_R_op = G_R_op + G_R_noise;

        % Birth noise
        max_meas = w_f;
        min_meas = -w_f;
        B_R_noise = min_meas + (max_meas-min_meas) .* rand(1,1);
        B_R_op = B_R_op + B_R_noise;

        % Mortality noise
        max_meas = w_m;
        min_meas = -w_m;
        M_R_noise = min_meas + (max_meas-min_meas) .* rand(1,1);
        M_R_op = M_R_op + M_R_noise;

        %Update the matrix A (synthetic data)
        A_cont =compute_A_continous(Pest_stages);
        sysc = ss(A_cont,[],eye(8),[]);
        sysd = c2d(sysc,1,'zoh');
        A_dis = sysd.A;

        %Update the matrix A(estimation data)
        A_cont_op =compute_A_continous(Pest_stages_op);
        sysc = ss(A_cont_op,[],eye(8),[]);
        sysd = c2d(sysc,1,'zoh');
        A_dis_op = sysd.A;
        
       L = [-x_est(1) -x_est(1) x_est(8) zeros(1,12) ;... %Egg (1)
        x_est(1) 0 0 -x_est(2) -x_est(2) zeros(1,10);... %L1 (2)
        zeros(1,3) x_est(2) 0 -x_est(3) -x_est(3) zeros(1,8);... %L2 (3)
        zeros(1,5) x_est(3) 0 -x_est(4) -x_est(4) zeros(1,6);... %L3 (4)
        zeros(1,7) x_est(4) 0 -x_est(5) -x_est(5) zeros(1,4);...%P (5)
        zeros(1,9) 0.5*x_est(5) 0 -x_est(6) zeros(1,3);...%AM (6)
        zeros(1,9) 0.5*x_est(5) 0 0  -x_est(7) -x_est(7) 0;%NMF (7)
        zeros(1,12) x_est(7) 0 -x_est(8)]; %MF (8)
    end
    
    %Store the variables
    ws_name=strcat(pwd,'\iterations\','ws',num2str(it));
    save(ws_name, 'x_hist', 'x_open_hist', 'x_est_hist')

end

%% PROCESS DATA

%Initialize indicators
RMSE_open = [];
RMSE_ekf = [];

r2_open = [];
r2_ekf = [];

%Compute values for all iterations
for ind=1:Number_iterations
    ws_name=strcat(pwd,'\iterations\','ws',num2str(ind));
    load(ws_name)
    
    temp_open = 0 ;
    temp_ekf = 0;
   
    temp2_open = 0 ;
    temp2_ekf = 0;
  
    temp3_open = 0 ;
    temp3_ekf = 0;
 
    
    for j = 1:Simulation_time
        
        %RMSE
        temp_open = temp_open + (x_hist(6,j)-x_open_hist(6,j))^2;
        temp_ekf = temp_ekf + (x_hist(6,j)-x_est_hist(6,j))^2; 
        
        %Chi squared fit
        temp2_open = temp2_open + ((x_hist(6,j)-x_open_hist(6,j))^2)/x_hist(6,j);
        temp2_ekf = temp2_ekf + ((x_hist(6,j)-x_est_hist(6,j))^2)/x_hist(6,j);

       
        %R2 fit
        temp3_open = temp3_open + (mean(x_open_hist(6,:))-x_open_hist(6,j))^2 ;
        temp3_ekf = temp3_ekf + (mean(x_est_hist(6,:))-x_est_hist(6,j))^2 ;
       
        
    end
    RMSE_open = [RMSE_open, sqrt(temp_open/Simulation_time)];
    RMSE_ekf = [RMSE_ekf, sqrt(temp_ekf/Simulation_time)];
  
    r2_open = [r2_open, 1 - (temp_open/temp3_open)];
    r2_ekf = [r2_ekf, 1 - (temp_ekf/temp3_ekf)];
   
end

%Average RMSE
Avg_open =mean(RMSE_open);
Avg_ekf = mean(RMSE_ekf);


%Average R2
Avg_r2_open = mean(r2_open);
Avg_r2_ekf = mean(r2_ekf);


%% PLOTS
t1 = datetime(2020,3,15,12,0,0);
t=t1+days(0:Simulation_time-1);

tt = 1:Simulation_time;

% Plot adult males population  

figure
plot(t,x_hist(6,:),'LineWidth',2);
hold on
plot(t,x_open_hist(6,:),'LineWidth',2);
hold on
plot(t,x_est_hist(6,:),'LineWidth',2);
legend("Adult males","Open loop adult males","EKF ",'Fontsize',15);
xlabel('Time [days]','Fontsize',15);
ylabel('Number of adult males','Fontsize',15);
title('Synthetic data PANTHEON 2020','Fontsize',20);
set(gca,'FontSize',15)

%% PLOTS
figure
boxplot([RMSE_open',RMSE_ekf'],'Labels',{'Open loop','EKF'})
title('RMSE')

figure
boxplot([r2_open',r2_ekf'],'Labels',{'Open loop','EKF'})
title('R2 fit')
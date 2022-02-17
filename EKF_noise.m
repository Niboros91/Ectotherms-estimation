function xhatOut = EKF_noise(iteration,A,L,C,N,x_ini,meas,measuring)
% Define storage for the variables that need to persist
% between time periods.
persistent xhat Q R P 

if iteration==1 %In the first iteration we create the 

    xhat=x_ini;
    P = eye(N)*1;
    R = eye(3)*2;     % Checked
    Q = eye(15)*0.005;
   
end

% Estimation Kalman filter
xhat = A*xhat;
P = A*P*transpose(A) + L*Q*transpose(L);
P=(P+P.')/2; %To force symmetry
P_pre=P;


if measuring == 1
    % Update Kalman Filter
    
    K =( P*C')/((C*P*C' + R));
    xhat=xhat+K*(meas-C*xhat);
    P=P-K*C*P;
end

% Post the results
xhatOut = xhat;
end
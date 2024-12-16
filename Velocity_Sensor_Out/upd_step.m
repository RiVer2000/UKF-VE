function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
%% Parameter Definition
%z_t - is the sensor data at the time step
%covarEst - estimated covar of the  state
%uEst - estimated mean of the state
%% %%%%%%%%%%%%%%%%%% My Implementation starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UKF Parameters
alpha = 0.001; 
beta = 2; 
kappa = 1; 


sigma_points = [];                      % Placeholder matrix for the sigma points 
sigma_points_prop = [];                 % Placeholder matrix for the propagated sigma pounts
n = 15;                                 % Number of Sigma Points
lambda = ((alpha^2)*(n + kappa)) - n;    
z_ut = [];                              % Placeholder matrix for updated z_ut
C_t = [];                               % Placeholder matrix for Cross Covariance
S_t = [];                               % Placeholder matrix for Covariance of the output

% Frame Transformation and Rotations
R_c_b = [0.707 -0.707  0; 
        -0.707 -0.707  0; 
             0      0 -1];

T = [-0.04; 0.0; -0.03];  
T_b_c =   inv([R_c_b T; 0 0 0 1]);  

R_bc = T_b_c(1:3,1:3); 

skew_m = [      0    0.0300   -0.0283;
          -0.0300         0   -0.0283;
           0.0283    0.0283         0];

% Diagonal matrix of noises
R_n = diag([0.06; 0.01; 0.04]);
% R_n = diag([0.04; 0.1; 0.04]);

% Cholesky Decomposition
covarChol = chol(covarEst, "lower");  

%Computing the sigma points
for i = 1 : 2*n + 1 
    if i == 1
        sigma_points(1 : 15, i) = uEst;  
    elseif i <= 16
        sigma_points(1 : 15, i) = uEst + (sqrt(n + lambda) * covarChol(:, i - 1)); 
    elseif i > 16 
        sigma_points(1 : 15, i) = uEst - (sqrt(n + lambda) * covarChol(:, rem(i, 16))); 
    end
end

%Propogate the sigma points
for i = 1 : 2*n +1 
    %Rotation matrix
    R = [cos(sigma_points(5, i))*cos(sigma_points(6, i))         cos(sigma_points(6, i))*sin(sigma_points(4, i))*sin(sigma_points(5, i)) - cos(sigma_points(4, i))*sin(sigma_points(6, i))          sin(sigma_points(4, i))*sin(sigma_points(6, i)) + cos(sigma_points(4, i))*cos(sigma_points(6, i))*sin(sigma_points(5, i));
         cos(sigma_points(5, i))*sin(sigma_points(6, i))         cos(sigma_points(4, i))*cos(sigma_points(6, i)) + sin(sigma_points(4, i))*sin(sigma_points(5, i))*sin(sigma_points(6, i))          cos(sigma_points(4, i))*sin(sigma_points(5, i))*sin(sigma_points(6, i)) - cos(sigma_points(6, i))*sin(sigma_points(4, i));
                               - sin(sigma_points(5, i))                                                                                   cos(sigma_points(5, i))*sin(sigma_points(4, i))                                                                                   cos(sigma_points(4, i))*cos(sigma_points(5, i))]; 

    %Computing linear velocity
    sigma_points_prop(1 : 3, i) = (R_c_b * transpose(R) * sigma_points(7 : 9, i)) - (R_c_b * skew_m * transpose(R_c_b) * z_t(4 : 6, 1));
end

%Computing the  z_ut
for i = 1 : 31 
    if (i == 1) 
        W_m = lambda / (n+lambda); 
        z_ut = W_m * sigma_points_prop(:, i); 
    else
        W_m = 1 / (2 * (n + lambda)); 
        z_ut = z_ut + (W_m * sigma_points_prop(:, i)); 
    end
end

%Computing the Cross Covariance
for i = 1 : 2*n + 1 
    if (i == 1) 
        W_c = (lambda / (n + lambda)) + (1 - alpha^2 + beta); 
        C_t = W_c * (sigma_points(1 : 15, i) - uEst) * (transpose(sigma_points_prop(:, i) - z_ut)); 
    else
        W_c = 1/(2 * (n + lambda)); 
        C_t = C_t + (W_c * (sigma_points(1 : 15, i) - uEst) * (transpose(sigma_points_prop(:, i) - z_ut))); 
    end
end

%Computing the Coavriance of the output
for i = 1 : 2*n + 1  
    if (i == 1) 
        W_c = (lambda / (n + lambda)) + (1 - alpha^2 + beta); 
        S_t = (W_c * (sigma_points_prop(:, i) - z_ut) * (transpose(sigma_points_prop(:, i) - z_ut)));
    else
        W_c = 1/(2 * (n + lambda)); 
        S_t = (S_t + (W_c * (sigma_points_prop(:, i) - z_ut) * (transpose(sigma_points_prop(:, i) - z_ut)))); 
    end
end

%Recomputing S_t by adding noise 
S_t = S_t + R_n;

%Computing the Filter Gain
K = C_t / S_t; 

%Computing the Filtered Mean and Covariance
uCurr = uEst + (K*(z_t(1 : 3, 1) - z_ut));
covar_curr = covarEst - (K * S_t * transpose(K)); 

%% %%%%%%%%%%%%% END %%%%%%%%%%%%%%%%    
end
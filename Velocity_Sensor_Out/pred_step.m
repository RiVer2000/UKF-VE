function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
%% Parameter Definition
% uPrev - is the mean of the prev state
%covarPrev - covar of the prev state
%angVel - angular velocity input at the time step
%acc - acceleration at the timestep
%dt - difference in time 

%% %%%%%%%%%%%%%%%% My Implementation starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State Constants
g = [0 ; 0; -9.81];
ng = ones(3,1)*0.01;           % Noise Gyroscope (From Tuning)
na = ones(3,1)*0.01;           % Noise Accelerometer (From Tuning)
nbg = ones(3,1)*0.01;          % Bias Gyroscope (From Tuning)
nba = ones(3,1)*0.01;          % Bias Accelerometer (From Tuning)
omega_m = angVel(1:3);          % Angular Velocity
a_m = acc(1:3);                 % Linear Acceleration


% UKF Parameters
alpha = 0.001;
kappa = 1;
beta = 2;

sigma_points = [];              % Placeholder matrix for the sigma points 
sigma_points_prop = [];         % Sigma Points Propogation matrix
n = 15;                         % State Dimension
n_q = 12;                       % Noise Dimension
n_prime = n + n_q;

lambda_prime = (alpha ^ 2) * (n + kappa) - n;


%% To Do
% create augmented state x_aug = [x ; q] with mean u = [u; 0] 
% covar_aug = [covar 0; 0 Q] cholsky decomposition


% Compute Sigma Points
uAug = [uPrev; zeros(12,1)];
covarAug = [covarPrev    zeros(15,12);
            zeros(12,15) diag([ng; na; nbg; nba])];



% Cholesky Decomposition 
covarChol = chol(covarAug, "lower");

% Compute Sigma Points
% The points are set such that the first column has the first sigma point 
% and the positive sigma points are in the 2nd to 28th column and the 
% negative are in the rest columns.
for i = 1: 2 * n_prime + 1
    if i==1
        sigma_points(1 : 27, i) = uAug;
    elseif i <=28
        sigma_points(1 : 27, i) = uAug + sqrt(n_prime + lambda_prime) * covarChol(:, i - 1);
    elseif i > 28
        sigma_points(1 : 27, i) = uAug - sqrt(n_prime + lambda_prime) * covarChol(:, rem(i, 28));
    end
end

%Propogating Sigma Points through the nonlinear function h
for i = 1 : 2 * n_prime + 1

    R = [cos(sigma_points(5, i))*cos(sigma_points(6, i))         cos(sigma_points(6, i))*sin(sigma_points(4, i))*sin(sigma_points(5, i)) - cos(sigma_points(4, i))*sin(sigma_points(6, i))          sin(sigma_points(4, i))*sin(sigma_points(6, i)) + cos(sigma_points(4, i))*cos(sigma_points(6, i))*sin(sigma_points(5, i));
         cos(sigma_points(5, i))*sin(sigma_points(6, i))         cos(sigma_points(4, i))*cos(sigma_points(6, i)) + sin(sigma_points(4, i))*sin(sigma_points(5, i))*sin(sigma_points(6, i))          cos(sigma_points(4, i))*sin(sigma_points(5, i))*sin(sigma_points(6, i)) - cos(sigma_points(6, i))*sin(sigma_points(4, i));
                               - sin(sigma_points(5, i))                                                                                   cos(sigma_points(5, i))*sin(sigma_points(4, i))                                                                                   cos(sigma_points(4, i))*cos(sigma_points(5, i))]; 

    G_inv = pinv([cos(sigma_points(5, i))*cos(sigma_points(6, i))         -sin(sigma_points(6, i))            0;
                  cos(sigma_points(5, i))*sin(sigma_points(6, i))          cos(sigma_points(6, i))            0;
                                         -sin(sigma_points(5, i))                               0            1])* R;

    sigma_points_prop(1 : 15, i) = sigma_points(1 : 15, i) + (dt*[sigma_points(7 : 9, i); ...
                                                        G_inv * (omega_m - sigma_points(10 : 12, i) - sigma_points(16 : 18, i)); ...
                                                        g + (R * (a_m - sigma_points(13 : 15, i) - sigma_points(19 : 21, i))); ...
                                                        sigma_points(22 : 24, i); ...
                                                        sigma_points(25 : 27, i)]); %Propogation of sigma points
end

% Compute the Predicted Mean 
for i = 1 : 55 
    if (i == 1) 
        W_m = lambda_prime / (n_prime + lambda_prime);     % Wm for the first sigma point
        uEst = W_m * sigma_points_prop(:, i);             % Mean for the first sigma point
    else
        W_m = 1 / (2 * (n_prime + lambda_prime));                % Wm for rest of the sigma points
        uEst = uEst + (W_m * sigma_points_prop(:, i));   % Mean for rest of the sigma points
    end
end


% Compute Covariance
for i = 1 : 55 %for loop for all the propgated sigma points
    if (i == 1) 
        W_c = (lambda_prime / (n_prime + lambda_prime)) + (1 - alpha^2 + beta); %Computing Wc for the first sigma point
        covarEst = W_c * (sigma_points_prop(:, i) - uEst) * transpose((sigma_points_prop(:, i) - uEst)); %Computing the Covariance for the first sigma point
    else 
        W_c = 1/(2 * (n_prime + lambda_prime)); %Computing Wc for rest of the sigma points
        covarEst = covarEst + (W_c * (sigma_points_prop(:, i) - uEst) * transpose((sigma_points_prop(:, i) - uEst))); %Computing the Covariance for rest of the sigma points
    end
end


%% %%%%%%%%%%%%% END %%%%%%%%%%%%%%%%
end


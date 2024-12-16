function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state

    %% To do C R Kalman gain Compute mean and covar
    % Using pose estimation from the optical flow
    Ct = [eye(6), zeros(6, 9)];

    Rt = eye(6) * 0.003;

    %Kalman gain 
    Kt = (covarEst * transpose(Ct))*pinv((((Ct * covarEst * transpose(Ct)) + Rt))); 

    % Current Mean
    uCurr = double(uEst + (Kt * (z_t - (Ct * uEst))));

    % Current Covariance
    covar_curr = double(covarEst - (Kt * Ct * covarEst));
    
end


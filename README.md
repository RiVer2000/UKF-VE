# State Estimation and Optical Flow based Tracking of a quadrotor.
## Introduction
In this project we develop an Unscented Kalman Filter (UKF) for our IMU driven quadcopter model to fuse the inertial data and the vision-based pose and velocity estimation. The UKF capturtes the non linearity better than Extended Kalman Filter (EKF) although is is computationally heavier. Unlike the EKF which uses Jacobian to deal with the non linearity, the UKF deals with nonlinearity by using a deterministic sampling approach to capture the mean and covariance estimates.
## Pose Estimation
We implement a vision based 3-D pose estimator which estimates position and orientation of the quadrotor based on AprilTags. The quadcopter was flown over a AprilTag mat of 108 tags set in 12x9 matrix. For each AprilTag we use 4 corners and 1 center point coordinates. We can also use the relationship between the coordinates in the world frame and the camera frame and for each set of corresponding points in world and camera frame. We have 5 coordinates for each Apriltag that we stack vertically and get the A matrix of size 10x9. Then in order to get Rotation and translation, we multiply the Homography matrix with inverse of the camera calibration matrix. We also need to normalize the Rotation and Translation to satify the constraints. Next we take the Singular Value Decomposition and from the obtained matrices we can obtain the Rotation Matrix. 
## Velocity Estimation
The features of two images are extracted which are used in the equations to calculate A(p) and B(p) which uses depth(z) values. Then the two matrixes A(p) and B(p) are combined together in one H matrix whose pseudo-inverse is multiplied with pdot which is multiplied with the adjoint matrix to get the estimated velocities. We also filter the values using Savitzky-Golay filtering to get cleaner graphs.
## System Model
The quadcopter’s state-space model includes:
- Position, Orientation, Linear Velocity: Key motion variables.
- Gyroscope and Accelerometer Bias: To account for sensor drift.
- Non-Linear Dynamics: Handled using deterministic sigma-point propagation.
### Key features of the model:
State-Space Augmentation: Incorporates IMU and vision data for fusion.
Non-Linear Orientation Representation: Utilizes Z-Y-X Euler angles for SO(3) transformations.
## Process Model:
- Gyroscope and accelerometer provide noisy estimates of angular velocity and linear acceleration. 
- Biases are modeled as Gaussian white noise.
## Update Model:
- Vision-based pose and velocity measurements update the UKF.
- Linear additive models are used for pose updates.

Results:
![image](https://github.com/user-attachments/assets/81622ec9-6e97-4ccf-ac3f-9a9eeb478f79)



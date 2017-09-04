## Unscented Kalman Filter (UKF)
- Code ompiles without errors with cmake and make.
- As tabulated below, px, py, vx, vy output coordinates have an RMSE that is less than the __threshold__: [.09, .10, 0.40, 0.30] 

  | Dataset | px        | py        | vx       | vy       |
  |---------|-----------|-----------|----------|----------|
  | 1       | 0.0734922 | 0.0819829 | 0.323838 | 0.184022 |
  | 2       | 0.0714487 | 0.0707447 | 0.253624 | 0.201682 |

- I tuned the process noise standard deviation of longitudinal acceleration, std_a, and the yaw acceleration, std_yawdd, to the
following values.

  | std_a (m/s^2) | std_yawdd  (rad/s^2) |
  |---------------|----------------------|
  | 1.5           | M_PI / 8             |

- I tuned the covariance matrix, P, by setting the elements of major diagonal to scaled down values of variances observed from the
datasets (for example, var_px = var_px / 100).

- I plotted the NIS (Normalized Innovation Sqaured) metric collected for both lidar and radar. The NIS values are greater than
5.991 about 3.81% of all lidar cases. For radar, the NIS values are greater than 7.815 about 4% of all cases.

- RMSE velocity is lower for UKF compare to that of EKF. This shows that UKF handles non-linear equations better than EKF. And, the use of
CTRV (Constant Turn Rate and Velocity) model enables more precise calculations.

- When turning off either radar data or lidar data, RMSE is greater than the threshold. Lidar provides more accurate readings 
compare to radar. The fact that lidar has higher resolution than radar explains the reason. On the other hand, radar data includes
the speed of moving object measured directly which lidar is not capable of measuring. Fusing both lidar and radar data and using UKF 
techniques to handle non-linear equations to predict and update the state of object resuts in improved results.

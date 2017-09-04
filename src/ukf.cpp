#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Set state dimension
  n_x_ = 5;

  // Set augmented state dimension
  n_aug_ = 7;

  // Define sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Set weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(1.0 / (2.0 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.setZero();

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_.setIdentity();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 8;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Instantaneous NIS for radar
  NIS_radar_ = 0.0;

  // Instantaneous NIS for laser
  NIS_laser_ = 0.0;

  // Set measurement dimension, lidar can measure px and py
  n_z_laser_ = 2;

  // Set measurement dimension, radar can measure r, phi, and r_dot
  n_z_radar_ = 3;

  // Measurement covariance matrix - laser
  R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  R_laser_ << (std_laspx_ * std_laspx_), 0,
              0, (std_laspy_ * std_laspy_);

  // Measurement covariance matrix - radar
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_ << (std_radr_ * std_radr_), 0, 0,
              0, (std_radphi_ * std_radphi_), 0,
              0, 0, (std_radrd_ * std_radrd_);
}

UKF::~UKF() {}

// -----------------------------------------------------------------------------
// @param {MeasurementPackage} meas_package: The latest measurement data of
// either radar or laser.
// -----------------------------------------------------------------------------
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) {
    // Initial state covariance matrix P
    // Obtained by observing standard deviation of given data sets
    P_ << 2.3075, 0, 0, 0, 0,
          0, 0.8298, 0, 0, 0,
          0, 0, 0.1765, 0, 0,
          0, 0, 0, 0.0239, 0,
          0, 0, 0, 0, 0.0015;

    // Initialize state.
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      double ro = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      // Radar velocity, rho_dot, is measured from the autonomous vehicle's perspective.
      // CTRV velocity, v, is tangential to the circle along which the bicycle travels.
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v = sqrt(vx * vx + vy * vy);
      v = (rho_dot < 0.0) ? -v : v;

      x_ << ro * cos(phi), ro * sin(phi), v, 0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    // Time when the state is true, in us
    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  // Compute the time elapsed between the current and previous measurements
  // Expressed in seconds
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  if (dt == 0.0) return;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  while (dt > 0.1) {
    const double dt_incr = 0.05;
    Prediction(dt_incr);
    dt -= dt_incr;
  }
  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  // Update the state and covariance matrices with the most recent raw measurements
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } else {
    // Laser updates
    UpdateLidar(meas_package);
  }

  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
}

// -----------------------------------------------------------------------------
// Predicts sigma points, the state, and the state covariance matrix.
// @param {double} delta_t the change in time (in seconds) between the last
// measurement and this one.
// -----------------------------------------------------------------------------
void UKF::Prediction(double delta_t) {
  // Constants
  const int n_sigma = 2 * n_aug_ + 1;    // 2n + 1 sigma points

  // Create augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // Create augmented covariance matrix, P_aug
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;

  // Add Q, the uncertainty caused by process noise, to P_aug
  P_aug(n_aug_ - 2, n_aug_ - 2) =  std_a_ * std_a_;
  P_aug(n_aug_ - 1, n_aug_ - 1) =  std_yawdd_ * std_yawdd_;

  // Calculate square root of P_aug
  MatrixXd A = P_aug.llt().matrixL();

  // Create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma);
  MatrixXd var = sqrt(lambda_ + n_aug_) * A;    // variance
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + var.col(i);
    Xsig_aug.col(i + n_aug_ + 1) = x_aug - var.col(i);
  }

  // Predict sigma points
  for (int i = 0; i < n_sigma; i++) {
    VectorXd x = Xsig_aug.col(i);
      
    // Add the longitudinal acceleration and yaw acceleration process noises
    // Calculation is same regardless of the yaw rate
    VectorXd v(n_x_);
    v <<  0.5 * delta_t * delta_t * cos(x(3)) * x(5),
          0.5 * delta_t * delta_t * sin(x(3)) * x(5),
          delta_t * x(5),
          0.5 * delta_t * delta_t * x(6),
          delta_t * x(6); 
    
    // Add the deterministic part of the CTRV process model
    // Avoid division by zero as well
    if (fabs(x(4)) >= 0.0001) {
      // Non-zero yaw rate - object turning
      v(0) += x(2) / x(4) * (sin(x(3) + x(4) * delta_t) - sin(x(3)));
      v(1) += x(2) / x(4) * (-cos(x(3) + x(4) * delta_t) + cos(x(3)));
      v(3) += x(4) * delta_t;
    } else {
      // Yaw rate is zero - object moving in straight line
      v(0) += x(2) * cos(x(3)) * delta_t;
      v(1) += x(2) * sin(x(3)) * delta_t;
    }
      
    // Predicted sigma points x_k+1 = x_k + v
    Xsig_pred_.col(i) = x.head(n_x_) + v;
  }

  // Predict state mean (vectorizing)
  x_ = Xsig_pred_ * weights_;
  
  //Predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sigma; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Normalize yaw angle to [-pi, pi)
    Tools::NormalizeAngle1(x_diff(3)); 

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

// -----------------------------------------------------------------------------
// Updates the state and the state covariance matrix using a laser measurement.
// @param {MeasurementPackage} meas_package
// -----------------------------------------------------------------------------
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // Constants
  const int n_sigma = 2 * n_aug_ + 1;    // 2n + 1 sigma points

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_laser_, n_sigma);

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_laser_);
  
  // Transform sigma points into measurement space
  z_pred.fill(0.0);
  for (int i = 0; i < n_sigma; i++) {
    double px = Xsig_pred_.col(i)(0);
    double py = Xsig_pred_.col(i)(1);
      
    // Create matrix for sigma points in measurement space
    Zsig.col(i)(0) = px;
    Zsig.col(i)(1) = py;
      
    // Calculate mean predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  // Predict measurement covariance matrix S
  // S_k+1 = weights * z_diff * z_diff_T + R
  MatrixXd S = MatrixXd(n_z_laser_, n_z_laser_);
  S.fill(0.0);
  for (int i = 0; i < n_sigma; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  // Add measurement noise covariance matrix, R, to S
  S += R_laser_;

  // Calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z_laser_);
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma; i++) {
    // State vector diff
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // Normalize yaw angle to [-pi, pi)
    Tools::NormalizeAngle1(x_diff(3)); 
      
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
      
    // Cross Correlation
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  // Calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  // Update state mean and covariance matrix
  // Incoming lidar measurements
  VectorXd z = meas_package.raw_measurements_;
  // Residual - diff in position between measurement and prediction: 2x1
  VectorXd z_diff = z - z_pred;
  // New estimate
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // NIS - Normalized Innovation Squared
  // Lidar has 2 degree of freedom => NIS range: 0.103 - 5.991 for 95% of cases
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  //cout << "NIS_laser " << NIS_laser_ << endl;
}

// -----------------------------------------------------------------------------
// Updates the state and the state covariance matrix using a radar measurement.
// @param {MeasurementPackage} meas_package
// -----------------------------------------------------------------------------
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // Constants
  const int n_sigma = 2 * n_aug_ + 1;    // 2n + 1 sigma points

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar_, n_sigma);

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_);
  
  // Transform sigma points into measurement space
  z_pred.fill(0.0);
  for (int i = 0; i < n_sigma; i++) {
    double px = Xsig_pred_.col(i)(0);
    double py = Xsig_pred_.col(i)(1);
    double v = Xsig_pred_.col(i)(2);
    double yaw = Xsig_pred_.col(i)(3);
    double rho = sqrt(px * px + py * py);
    double phi = (px == 0.0) ? 0.0 : atan2(py, px);
    double rho_dot = 0.0;
    if (fabs(rho) >= 0.0001) {
      rho_dot = (px * cos(yaw) * v + py * sin(yaw) * v) / rho;
    }
      
    // Create matrix for sigma points in measurement space
    Zsig.col(i)(0) = rho;
    Zsig.col(i)(1) = phi;
    Zsig.col(i)(2) = rho_dot;
      
    // Calculate mean predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  // Predict measurement covariance matrix S
  // S_k+1 = weights * z_diff * z_diff_T + R
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);
  S.fill(0.0);
  for (int i = 0; i < n_sigma; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // Normalize phi angle to [-pi, pi)
    Tools::NormalizeAngle1(z_diff(1)); 
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  // Add measurement noise covariance matrix, R, to S
  S += R_radar_;

  // Calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma; i++) {
    // State vector diff
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // Normalize yaw angle to [-pi, pi)
    Tools::NormalizeAngle1(x_diff(3)); 
      
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // Normalize phi angle to [-pi, pi)
    Tools::NormalizeAngle1(z_diff(1)); 
      
    // Cross Correlation
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  // Calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  // Update state mean and covariance matrix
  // Incoming radar measurements
  VectorXd z = meas_package.raw_measurements_;
  // Residual - diff in position between measurement and prediction: 3x1
  VectorXd z_diff = z - z_pred;
  // Normalize phi angle to [-pi, pi)
  Tools::NormalizeAngle1(z_diff(1)); 
  // New estimate
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // NIS - Normalized Innovation Squared
  // Radar has 3 degree of freedom => NIS range: 0.35 - 7.81 for 95% of cases
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  //cout << "NIS_radar " << NIS_radar_ << endl;
}

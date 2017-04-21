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

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_ << 0, 0, 0, 0, 0;

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ <<
    1, 0, 0, 0, 0,
    0, 1, 0, 0, 0,
    0, 0, 1000, 0, 0,
    0, 0, 0, 1000, 0,
    0, 0, 0, 0, 1000;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.7;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.25;

  YAW_ACCEL_MAX_ = 0.8;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_; 

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  // time when the state is true, in us
  time_us_ = 0.0;

  // Sigma point spreading parameter
  float lambda_ = 3 - n_x_;

  // Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.0);
  // set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  // the current NIS for radar
  NIS_radar_ = 0.0;

  // the current NIS for laser
  NIS_laser_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if ( !is_initialized_ ) {
    // First measurement
    cout << "UKF: " << endl;

    if ( meas_package.sensor_type_ == MeasurementPackage::RADAR ) {
      //Convert radar from polar to cartesian coordinates and initialize state.
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      // can you guess small values for v, phi and phi_dot?
      x_ << px, py, 0, 0, 0;
    } else if ( meas_package.sensor_type_ == MeasurementPackage::LASER ) {
      //Initialize state.
      x_ << meas_package.raw_measurements_[0],
	    meas_package.raw_measurements_[1],
            0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;

    // Done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  if ( meas_package.sensor_type_ == MeasurementPackage::RADAR  && !use_radar_) {
    return void();
  }
  if ( meas_package.sensor_type_ == MeasurementPackage::LASER  && !use_laser_) {
    return void();
  }
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  // Compute the time elapsed between the current and previous measurements
  // delta_t is in seconds
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  // normalize angle
  //x_(3) = tools_.NormalizeAngle(x_(3));
  // limit angle acceleration
  //if (x_(4)> 0) x_(4) = min(x_(4), YAW_ACCEL_MAX_);
  //if (x_(4)< 0) x_(4) = max(x_(4), -YAW_ACCEL_MAX_);
  // limit max speed
  //if (x_(2)> 0)  x_(2) = min(x_(2), 6.0);
  //if (x_(2)< 0)  x_(2) = max(x_(2), -6.0);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if ( meas_package.sensor_type_ == MeasurementPackage::RADAR ) {
    // Radar updates
    UpdateRadar(meas_package);
  } else {
    // Laser updates
    UpdateLidar(meas_package);
  }

  // normalize angle
  //x_(3) = tools_.NormalizeAngle(x_(3));
  // limit angle acceleration
  //if (x_(4)> 0) x_(4) = min(x_(4), YAW_ACCEL_MAX_);
  //if (x_(4)< 0) x_(4) = max(x_(4), -YAW_ACCEL_MAX_);
  // limit max speed
  //if (x_(2)> 0)  x_(2) = min(x_(2), 6.0);
  //if (x_(2)< 0)  x_(2) = max(x_(2), -6.0);

  // Print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++) {
      Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //create matrix with predicted sigma points as columns
  //Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  //predict sigma points
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yaw_p = tools_.NormalizeAngle(yaw_p);
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
  //predicted state mean
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    x_diff(3) = tools_.NormalizeAngle(x_diff(3));
    // limit angle acceleration
    if (x_diff(4)> 0) x_diff(4) = min(x_diff(4), YAW_ACCEL_MAX_);
    if (x_diff(4)< 0) x_diff(4) = max(x_diff(4), -YAW_ACCEL_MAX_);

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  x_ = x;
  P_ = P;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  VectorXd z = meas_package.raw_measurements_;
  MatrixXd H_ = MatrixXd(2, n_x_);
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();

  MatrixXd S = H_ * P_ * Ht + R_laser_;
  // Check if S is invertable, if not, don't update
  Eigen::FullPivLU<MatrixXd> lu(S);
  if (!lu.isInvertible()) {
    return void();
  }
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  NIS_laser_ = y.transpose() * Si * y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  //create matrix for sigma points in measurement space
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = tools_.NormalizeAngle(z_diff(1));
 
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  // TODO: move to init
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_radr_*std_radr_, 0,                       0,
       0,                   std_radphi_*std_radphi_, 0,
       0,                   0,                       std_radrd_*std_radrd_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = tools_.NormalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = tools_.NormalizeAngle(x_diff(3));
    // limit angle acceleration
    if (x_diff(4)> 0) x_diff(4) = min(x_diff(4), YAW_ACCEL_MAX_);
    if (x_diff(4)< 0) x_diff(4) = max(x_diff(4), -YAW_ACCEL_MAX_);

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Check if S is invertable, if not, don't update
  Eigen::FullPivLU<MatrixXd> lu(S);
  if (!lu.isInvertible()) {
    return void();
  }
  MatrixXd Si = S.inverse();
  // Kalman gain K;
  MatrixXd K = Tc * Si;

  // Residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  // Angle normalization
  z_diff(1) = tools_.NormalizeAngle(z_diff(1));

  // Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  NIS_radar_ = z_diff.transpose() * Si * z_diff;
}




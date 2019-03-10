#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  /* init H_laser 
   * Lession25, section11 and 12 */
  H_laser_ << 1, 0, 0, 0, 
              0, 1, 0, 0;
  

  /* init Hj_
   * lession25, section19 */
  Hj_ << 1, 1, 0, 0,
         1, 1, 0, 0,
         1, 1, 1, 1;
  
  /* Set the initial transition matrix F_
   * Lession25 section9 and 13*/
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  /* process uncertainty matrix
   * Lession25, section7   */
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
  
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
      

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    /* Lession25, secction15 */
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      /* Lession25, section#15, 02:17 of video */
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2]; 
      
      /* Lession25, section#3, 01:03 of video */
      ekf_.x_(0) = rho * cos(phi);        
      ekf_.x_(1) = rho * sin(phi);
      ekf_.x_(2) = rho_dot * cos(phi);
      ekf_.x_(3) = rho_dot * sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_(0) = measurement_pack.raw_measurements_[0];
      ekf_.x_(1) = measurement_pack.raw_measurements_[1];
      ekf_.x_(2) = 0;
      ekf_.x_(3) = 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update      
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float noise_ax = 9; 
  float noise_ay = 9;
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

  /* add dt to F matrix, lession25, secction9 */
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  /* update convariance matrix Q, lession25, section10 */
  float dt_pow2 = dt * dt;
  float dt_pow3 = dt_pow2 * dt;
  float dt_pow4 = dt_pow3 * dt;
    
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             1, 0, 1, 0,
             0, 1, 0, 1;
     
  ekf_.Q_(0, 0) = dt_pow4/4 * noise_ax;
  ekf_.Q_(0, 2) = dt_pow3/2 * noise_ax;
  ekf_.Q_(1, 1) = dt_pow4/4 * noise_ay;
  ekf_.Q_(1, 3) = dt_pow3/2 * noise_ay;
  ekf_.Q_(2, 0) = dt_pow3/2 * noise_ax;
  ekf_.Q_(2, 2) = dt_pow2   * noise_ax;
  ekf_.Q_(3, 1) = dt_pow3/2 * noise_ay;
  ekf_.Q_(3, 3) = dt_pow2   * noise_ay;
  
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    Tools tools;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

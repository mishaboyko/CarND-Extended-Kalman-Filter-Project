#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

Eigen::IOFormat CleanMatrix(1, 0, ", ", "\n", "[", "]");
Eigen::IOFormat CleanVector(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "; ", "", "", "", ";");


/**
 * Constructor
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // first measurement
  cout << "EKF: " << endl;

  P_init_ = MatrixXd(4, 4);
  P_init_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 10, 0,
            0, 0, 0, 10;

  F_init_laser_ = MatrixXd(4, 4);
  F_init_laser_ <<  1, 0, 1, 0,
                    0, 1, 0, 1,
                    0, 0, 1, 0,
                    0, 0, 0, 1;

  // measurement (noise) covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement (noise) covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // laser measurement matrix
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // acceleration noise
  noise_ax = 9;
  noise_ay = 9;
}

/**
 * Destructor
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /**
   * Initialization
   */
  if (!is_initialized_) {

    ekf_.Init(P_init_, F_init_laser_, H_laser_, R_laser_, R_radar_);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      float x, y, v_x, v_y;
      
      x = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
      y = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
      v_x = measurement_pack.raw_measurements_[2] * cos(measurement_pack.raw_measurements_[1]);
      v_y = measurement_pack.raw_measurements_[2] * sin(measurement_pack.raw_measurements_[1]);

      ekf_.x_ << x, y, v_x, v_y;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  // compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  previous_timestamp_ = measurement_pack.timestamp_;

  // Update the state transition matrix F according to the new elapsed time.
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Update the process noise covariance matrix.
  ekf_.Q_ <<  dt_4/4*noise_ax,  0,                dt_3/2*noise_ax,  0,
              0,                dt_4/4*noise_ay,  0,                dt_3/2*noise_ay,
              dt_3/2*noise_ax,  0,                dt_2*noise_ax,    0,
              0,                dt_3/2*noise_ay,  0,                dt_2*noise_ay;

  /**
   * Prediction
   */
  ekf_.Predict();

  /**
   * Update
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.UpdateRadar(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.UpdateLaser(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "\nx_:\n" << "(" << ekf_.x_.format(CleanVector) << " )" << endl;
  cout << "\nP_:\n" << ekf_.P_.format(CleanMatrix) << endl;
}

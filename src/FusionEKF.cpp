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
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // first measurement
  cout << "EKF: " << endl;

  x_init_ = VectorXd(4);
  x_init_ << 1, 1, 1, 1;

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

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement (noise) covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement (noise) covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  // measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  /**
  * TODO: MBoiko: Hj_ << jacobian Matrix;
  */

  // set the acceleration noise
  noise_ax = 9;
  noise_ay = 9;

  Q_init = MatrixXd(4, 4);
  Q_init << 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;
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

    ekf_.Init(x_init_, P_init_, F_init_laser_, H_laser_, R_laser_, R_radar_, Q_init);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      float x, y;
      x = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
      y = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);

      ekf_.x_ << x, y, 0, 0;


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

  /**
   * Prediction
   */

  // Update the state transition matrix F according to the new elapsed time.
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Update the process noise covariance matrix.
  ekf_.Q_ <<  dt_4/4*noise_ax,  0,                dt_3/2*noise_ax,  0,
              0,                dt_4/4*noise_ay,  0,                dt_3/2*noise_ay,
              dt_3/2*noise_ax,  0,                dt_2*noise_ax,    0,
              0,                dt_3/2*noise_ay,  0,                dt_2*noise_ay;

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
    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
      // Laser updates
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "\nx_:\n" << "(" << ekf_.x_.format(CleanVector) << " )" << endl;
  cout << "\nP_:\n" << ekf_.P_.format(CleanMatrix) << endl;
}

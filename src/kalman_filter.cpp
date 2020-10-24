#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_laser_in, MatrixXd &R_radar_in) {
  SetInitValues();

  P_ = P_in;
  F_ = F_in;
  H_laser = H_in;
  R_laser = R_laser_in;
  R_radar = R_radar_in;
}

void KalmanFilter::SetInitValues(){
  x_ = VectorXd(4);
  x_ << 0, 0, 0, 0;

  Q_ = MatrixXd(4, 4);
  Q_ << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;

}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::UpdateLaser(const VectorXd &z) {

  // z and y in cartesian coordinate system
  Eigen::VectorXd y = z - H_laser * x_; // where H_ * x_ == predicted z

  UpdateCommon(y, H_laser, R_laser);
}

void KalmanFilter::UpdateRadar(const VectorXd &z) {

  // z and y in polar coordinate system
  Eigen::VectorXd y = z - tools.CartesianToPolar(x_); // where h(x') == predicted z

  // Normalize theta angles difference
  if (y[1] > M_PI) {
    y[1] = y[1] - 2 * M_PI;
  }
  else if (y[1] < -M_PI){
    y[1] = y[1] + 2 * M_PI;
  }

  // A matrix to convert from polar to cartesian
  H_radar = tools.CalculateJacobian(x_);

  UpdateCommon(y, H_radar, R_radar);
}

void KalmanFilter::UpdateCommon(Eigen::VectorXd y, Eigen::MatrixXd H, Eigen::MatrixXd R){

  Eigen::MatrixXd S = H * P_ * H.transpose() + R;
  // Calculate Kalman gain
  Eigen::MatrixXd K = (P_ * H.transpose()) * S.inverse();
  // new estimate
  x_ = x_ + (K * y);
  MatrixXd ident_matrix = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (ident_matrix - K * H) * P_;
}
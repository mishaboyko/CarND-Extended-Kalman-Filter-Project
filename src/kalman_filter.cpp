#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_laser_in, MatrixXd &R_radar_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_laser = H_in;
  R_laser = R_laser_in;
  R_radar = R_radar_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

  Eigen::VectorXd y = z - H_laser * x_; // where H_ * x_ == predicted z (z_pred)
  Eigen::MatrixXd S = H_laser * P_ * H_laser.transpose() + R_laser;
  // Calculate Kalman gain
  Eigen::MatrixXd K = (P_ * H_laser.transpose()) * S.inverse();
  // new estimate
  x_ = x_ + (K * y);
  MatrixXd ident_matrix = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (ident_matrix - K * H_laser) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  // y = z - h(x')
  Eigen::VectorXd y = z - tools.CartesianToPolar(x_); // where h(x_) == predicted z (z_pred)

  // Normalize theta angles difference
  if (y(1) > M_PI) {
    y(1) = y(1) - 2 * M_PI;
  }
  else if (y(1) < -M_PI){
    y(1) = y(1) + 2 * M_PI;
  }

  H_radar = tools.CalculateJacobian(x_);

  Eigen::MatrixXd S = H_radar * P_ * H_radar.transpose() + R_radar;
  cout << "s\n" << S << endl;
  // Calculate Kalman gain
  Eigen::MatrixXd K = (P_ * H_radar.transpose()) * S.inverse();
  // new estimate
  x_ = x_ + (K * y);
  MatrixXd ident_matrix = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (ident_matrix - K * H_radar) * P_;


  /*
  These three steps (initialize, predict, update) plus calculating RMSE encapsulate the entire extended Kalman filter project.
  */
}

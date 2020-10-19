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
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_laser = F_in;
  H_laser = H_in;
  R_laser = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_laser * x_;
  P_ = F_laser * P_ * F_laser.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  Eigen::VectorXd y = z - H_laser * x_; // where H_ * x_ == predicted z
  Eigen::MatrixXd S = H_laser * P_ * H_laser.transpose() + R_laser;
  Eigen::MatrixXd K = (P_ * H_laser.transpose()) * S.inverse();

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd ident_matrix = MatrixXd::Identity(x_size, x_size);
  P_ = (ident_matrix - K * H_laser) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  /*
  These three steps (initialize, predict, update) plus calculating RMSE encapsulate the entire extended Kalman filter project.
  */
}

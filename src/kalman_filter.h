#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include "tools.h"

class KalmanFilter {
 public:
  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param P_in Initial state covariance
   * @param F_in Transition matrix
   * @param H_in Measurement matrix
   * @param R_in Measurement covariance matrix
   * @param Q_in Process covariance matrix
   */
  void Init(Eigen::MatrixXd &P_in, Eigen::MatrixXd &F_in,
            Eigen::MatrixXd &H_in, Eigen::MatrixXd &R_laser_in, Eigen::MatrixXd &R_radar_in);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   */
  void Predict();

  /**
   * Updates the state by using _standard_ Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateLaser(const Eigen::VectorXd &z);

  /**
   * Updates the state by using _Extended_ Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateRadar(const Eigen::VectorXd &z);

  // object mean state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transition matrix ( for both, laser and radar)
  Eigen::MatrixXd F_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

  private:

  // set initial dummy values. These will be overwritten every sensor measurement
  void SetInitValues();

  // Update variables, which are common for EKF and KF (a.k.a. Radar and Laser measurements)
  void UpdateCommon(Eigen::VectorXd y, Eigen::MatrixXd H, Eigen::MatrixXd R);

  Tools tools;

  // measurement matrix for the laser
  Eigen::MatrixXd H_laser;

  // measurement matrix for the radar
  Eigen::MatrixXd H_radar;

  // measurement covariance matrix for the laser noise
  Eigen::MatrixXd R_laser;

  // measurement covariance matrix for the radar noise
  Eigen::MatrixXd R_radar;
};

#endif // KALMAN_FILTER_H_

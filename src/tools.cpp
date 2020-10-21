#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
   VectorXd rmse(4);
   rmse << 0,0,0,0;

   // Validity checking of the input:
   //  * the estimation vector size should not be zero
   //  * the estimation vector size should equal ground truth vector size
   
   if (estimations.size() != ground_truth.size() || estimations.size() == 0 || ground_truth.size() == 0) {
      cout << "Invalid matrices of estimation or the ground truth data" << endl;
      return rmse;
   }

   // sum up all square residuals

   for(unsigned int i = 0; i < estimations.size(); ++i){
      VectorXd residual = estimations[i] - ground_truth[i];

      // coefficient-wise multiplication
      residual = residual.array()*residual.array();

      rmse += residual;
   }

   // calculate the mean
   rmse = rmse / estimations.size();

   // calculate the squared root
   rmse = rmse.array().sqrt();

   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Input: a 4-D Vector x_state(x, y, v_x, v_y)
   * Output: 3x4 Jacobian Measurement matrix
   */

  Eigen::MatrixXd H_jacobian(3,4);

   // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute a set of terms
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  // check division by zero
  if (fabs(c1) < 0.0001) {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    exit(1);
  }

  // compute the Jacobian matrix
  H_jacobian << (px/c2),               (py/c2),                0,       0,
               -(py/c1),               (px/c1),                0,       0,
               py*(vx*py - vy*px)/c3,  px*(px*vy - py*vx)/c3,  px/c2,   py/c2;

  return H_jacobian;
}

VectorXd Tools::CartesianToPolar(const VectorXd& x_state) {
   // Range (rho). Radial distance from origin.
   float ro;
   // Bearing. Angle between rho and x-axis
   float theta;
   // radial velocity. Change of rho (range rate)
   float ro_dot;

   // recover state parameters
   float px = x_state(0);
   float py = x_state(1);
   float vx = x_state(2);
   float vy = x_state(3);

   /** 
                                                                       
   TODO:
   One other important point when calculating $y$ with radar sensor data: 
   the second value in the polar coordinate vector is the angle $\phi$. 
   You'll need to make sure to normalize $\phi$ in the $y$ vector 
   so that its angle is between $-\pi$ and $\pi$; in other words, 
   add or subtract $2\pi$ from $\phi$ until it is between $-\pi$ and $\pi$.
   
   
   */

   ro = sqrt(px*px+py*py);
   theta = atan2(py, px);
   ro_dot = (px*vx + py*vy)/ro;

   Eigen::VectorXd polar_space = VectorXd(3);
   polar_space << ro, theta, ro_dot;

   return polar_space;
}

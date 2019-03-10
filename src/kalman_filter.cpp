#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

const double PI = 3.141592653589793238463;
const double EPS = 0.0001;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state, Lession25-section9, Lession22-Section24
   */
  x_ = F_ * x_;
  
  /* Lession25- section10, it is a little bit different with lession22-section24, add Q */
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_; 
  
}

void KalmanFilter::CommonUpdateForLidarAndRadar(const VectorXd& y) {
  MatrixXd Ht     = H_.transpose();
  MatrixXd S      = H_*P_*Ht + R_;
  MatrixXd Si     = S.inverse();
  MatrixXd PHt    = P_*Ht;
  MatrixXd K      = PHt*Si;

  //new estimate
  x_          = x_ + (K*y);
  long x_size = x_.size();
  MatrixXd I  = MatrixXd::Identity(x_size, x_size);
  P_          = (I-K*H_) * P_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd y = z - (H_ * x_); 
  
  CommonUpdateForLidarAndRadar(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
    
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);
    
  const double rho = sqrt(px*px + py*py);
  const double phi = fabs(px) > EPS ? atan2(py, px) : 0.0; // atan2 returns values between -pi and pi
  const double rho_dot = fabs(rho) > EPS ? (px*vx + py*vy) / rho : 0.0;

  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  VectorXd y = z-z_pred;

  /* correction for phi in y to be not in the range -pi to pi
   * Errored on this place to make RMSE wrong. */
  y(1) = y(1)>PI ? y(1)-2*PI
       : y(1)<-PI ? y(1)+2*PI
       : y(1);

  CommonUpdateForLidarAndRadar(y);
}

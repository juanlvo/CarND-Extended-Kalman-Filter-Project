#include "kalman_filter.h"
#include <cmath>
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

void KalmanFilter::Predict(){

	x_ = F_ * x_ ; // There is no external motion, so, we do not have to add "+u"
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	  VectorXd z_pred = H_ * x_;
	  VectorXd y = z - z_pred;
	  MatrixXd Ht = H_.transpose();
	  MatrixXd S = H_ * P_ * Ht + R_;
	  MatrixXd Si = S.inverse();
	  MatrixXd PHt = P_ * Ht;
	  MatrixXd K = PHt * Si;
	  //new estimate
	  x_ = x_ + (K * y);
	  long x_size = x_.size();
	  MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  // Recalculate x object state to rho, theta, rho_dot coordinates

	  float rho_pred    =  pow(pow(x_[0],2) + pow(x_[1],2),0.5);
	  float phi_pred    =  0.0001;
	  if (fabs(x_[0]) > 0.0001) {
	    phi_pred  = atan2(x_[1],x_[0]);
	  }
	  float rhodot_pred = 0.0001;
	  if (fabs(rho_pred) > 0.0001) {
	    rhodot_pred = (x_[0]*x_[2] + x_[1]*x_[3]) / rho_pred;
	  }
	  VectorXd z_pred(3);
	  z_pred << rho_pred, phi_pred, rhodot_pred;

	  VectorXd y = z - z_pred;

	  //normalizing phi
	  while (y[1] > M_PI) {
		  y[1] = y[1]-M_PI;
	  }
	  /*while (y[1] < -(M_PI)) {
		  y[1] = y[1]+M_PI;
	  }*/

	  MatrixXd Ht = H_.transpose();

	  MatrixXd S = H_ * P_ * Ht + R_;
	  MatrixXd Si = S.inverse();
	  MatrixXd PHt = P_ * Ht;
	  MatrixXd K = PHt * Si;

	  //new estimate
	  x_ = x_ + (K * y);

	  while (x_[1] > (8*M_PI)) {
		  x_[1] = x_[1]-(8*M_PI);
	  }


	  long x_size = x_.size();
	  MatrixXd I = MatrixXd::Identity(x_size, x_size);
	  P_ = (I - K * H_) * P_;
}

// Step is a function execute every time for update the Kalman Filter
void KalmanFilter::Step(const VectorXd &y) {
   MatrixXd Ht = H_.transpose();
   MatrixXd S = H_ * P_ * Ht + R_;
   MatrixXd Si = S.inverse();
   MatrixXd PHt = P_ * Ht;
   MatrixXd K = PHt * Si;

   //new estimate
   x_ = x_ + (K * y);
   long x_size = x_.size();
   MatrixXd I = MatrixXd::Identity(x_size, x_size);
   P_ = (I - K * H_) * P_;
}


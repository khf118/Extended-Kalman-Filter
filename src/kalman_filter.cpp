#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
    TODO:
    * predict the state
    */
    cout << "x_" << x_ << endl;
    cout << "F_" << F_ << endl;
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
    MatrixXd y_, S_, K_;
    y_ = z - (H_ * x_);
    S_ = (H_ * P_ * H_.transpose()) + R_;
    K_ = P_ * H_.transpose() * S_.inverse();
    x_ = x_ + (K_ * y_);
    P_ = (MatrixXd::Identity(4,4) - (K_ * H_)) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    MatrixXd y_, S_, K_, Hj, Fj;
    float px = x_(0);
    float py = x_(1);
    
    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    if(fabs(c1) > 0.0001){
        VectorXd xpolar_(3);
        if (fabs(x_(0)) < 0.0001) {
            cout << "Devision by zero." << endl;
            return;
        }
        double c= sqrt( pow(x_(0),2) + pow(x_(1),2));
        xpolar_ << c, atan2(x_(1),x_(0)), (( x_(0) * x_(2) + x_(1) * x_(3)) / c);
        y_ = z - xpolar_;
        y_(1) = fmod(y_(1) + 3.14,6.28) - 3.14;
        S_ = (H_ * P_ * H_.transpose()) + R_;
        
        if (S_.determinant() != 0) {
            K_ = P_ * H_.transpose() * S_.inverse();
            x_ = x_ + (K_ * y_);
            P_ = (MatrixXd::Identity(4,4) - (K_ * H_)) * P_;
        }
    }
}

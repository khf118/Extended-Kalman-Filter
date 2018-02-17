#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
  
  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
    0, 1, 0, 0;
  ekf_.H_ = H_laser_;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
      
    ekf_.x_ = VectorXd(4);

    previous_timestamp_ = measurement_pack.timestamp_;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        VectorXd z(4);
        z= measurement_pack.raw_measurements_;
        
        double px = z(0) * cos(z(1));
        double py = z(0) * sin(z(1));
        ekf_.x_ << px, py, 1, 1;
        ekf_.F_ = MatrixXd(4,4);
        ekf_.F_ << 1,0,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1;
        ekf_.P_ = MatrixXd(4,4);
        ekf_.P_ << 10000,0,0,0,
        0,10000,0,0,
        0,0,10000,0,
        0,0,0,10000;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
        double x = measurement_pack.raw_measurements_(0);
        double y = measurement_pack.raw_measurements_(1);
        
        ekf_.x_ << x, y, 1, 1;
        ekf_.F_ = MatrixXd(4,4);
        ekf_.F_ << 1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1;
        ekf_.P_ = MatrixXd(4,4);
        ekf_.P_ << 10000,0,0,0,
        0,10000,0,0,
        0,0,10000,0,
        0,0,0,10000;
    }
      
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
    //compute the time elapsed between the current and previous measurements
    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds
    
    previous_timestamp_ = measurement_pack.timestamp_;
    double dt_2 = dt * dt;
    double dt_3 = dt_2 * dt;
    double dt_4 = dt_3 * dt;
    
    ekf_.F_ << 1,0,dt,0,
    0,1,0,dt,
    0,0,1,0,
    0,0,0,1;
    
    int noise_ax = 9;
    int noise_ay = 9;
    
    //set the process covariance matrix Q
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
    0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
    dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
    0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
  ekf_.Predict();
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      
      // Radar updates
      ekf_.R_ = R_radar_;
      
      VectorXd z = measurement_pack.raw_measurements_;
      double px = z(0) * cos(z(1));
      double py = -1 * z(0) * sin(z(1));
      VectorXd z2(4);
      z2 << px,py,2,3;
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.H_ = Hj_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
      // Laser updates
      ekf_.R_ = R_laser_;
      ekf_.H_ = H_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }
}

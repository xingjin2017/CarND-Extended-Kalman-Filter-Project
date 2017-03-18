#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

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
  */
  //create a 4D state vector, we don't know yet the values of the x state
  ekf_.x_ = VectorXd(4);

  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1000, 0,
    0, 0, 0, 1000;

  //measurement covariance
  ekf_.R_laser_ = R_laser_;
  ekf_.R_radar_ = R_radar_;

  //measurement matrix
  ekf_.H_laser_ = MatrixXd(2, 4);
  ekf_.H_laser_ << 1, 0, 0, 0,
    0, 1, 0, 0;

  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
    0, 1, 0, 1,
    0, 0, 1, 0,
    0, 0, 0, 1;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

VectorXd FusionEKF::NormalizeMeasurements(const MeasurementPackage &measurement_pack) {
  static float pi = static_cast<float>(3.14159265358979323846);
  static float two_pi = 2*pi;

  if (measurement_pack.sensor_type_ != MeasurementPackage::RADAR) {
    return measurement_pack.raw_measurements_;
  }
  VectorXd normalized(3);
  float ro = measurement_pack.raw_measurements_[0];
  float phi = measurement_pack.raw_measurements_[1];
  float ro_dot = measurement_pack.raw_measurements_[2];

  // Normalize phi to between [-pi, pi]
  while (abs(phi) > pi) {
    if (phi < 0) {
      phi += two_pi;
    } else {
      phi -= two_pi;
    }
  }

  normalized << ro, phi, ro_dot;
  return normalized;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  VectorXd normalized_measurements = NormalizeMeasurements(measurement_pack);

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
    cout << "EKF: " << endl;
    // ekf_.x_ = VectorXd(4);
    // ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = normalized_measurements[0];
      float phi = normalized_measurements[1];
      float ro_dot = normalized_measurements[2];
      ekf_.x_ << ro*cos(phi), ro*sin(phi), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << normalized_measurements[0],
	normalized_measurements[1], 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;
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
     * Use sigma_ax = 9 and sigma_ay = 9 for your Q matrix.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR
      && ekf_.x_.squaredNorm() < 0.0001) {
    // For radar, to avoid Hj zero dividing case, just skip the measurement.
    return;
  }
  
  float sigma_ax = 9;
  float sigma_ay = 9;

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*sigma_ax, 0, dt_3/2*sigma_ax, 0,
    0, dt_4/4*sigma_ay, 0, dt_3/2*sigma_ay,
    dt_3/2*sigma_ax, 0, dt_2*sigma_ax, 0,
    0, dt_3/2*sigma_ay, 0, dt_2*sigma_ay;

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
    cout << "ekf_.UpdateEKF" << endl;
    ekf_.UpdateEKF(normalized_measurements);
  } else {
    // Laser updates
    cout << "ekf_.Update" << endl;
    ekf_.Update(normalized_measurements);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
* Initializes Unscented Kalman filter
*/
UKF::UKF() {

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);    

    // Fixed measurement noise values provided by the sensor manufacturer:

    std_laspx_ = 0.15; // Laser measurement noise standard deviation position1 in m
    std_laspy_ = 0.15; // Laser measurement noise standard deviation position2 in m
    std_radr_ = 0.3; // Radar measurement noise standard deviation radius in m
    std_radphi_ = 0.03; // Radar measurement noise standard deviation angle in rad
    std_radrd_ = 0.3; // Radar measurement noise standard deviation radius change in m/s

    /**
    TODO:Complete the initialization. See ukf.h for other member properties.
    */
    
    std_a_ = 0.2; // Process noise standard deviation longitudinal acceleration in m/s^2

    std_yawdd_ = M_PI/8; // Process noise standard deviation yaw acceleration in rad/s^2

    is_initialized_ = false; // Used for init call to ProcessMeasurment

    time_us_ = 0; // Previous timestep
      
    n_x_ = 5; // State dimension

    n_aug_ = 7; // Augmented state dimension

    lambda_ = 3- n_x_; //Sigma point spreading parameter

    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ +1); // Init vector for sigma point predictions

    weights_ = VectorXd(2 * n_aug_ + 1); // Init vector for weights


    NIS_radar_ = 0.; //init Normalized Innovation Squared(NIS) for radar

    NIS_laser_ = 0.; //init Normalized Innovation Squared(NIS) for laser

}

UKF::~UKF() {}

/**
* @param {MeasurementPackage} meas_package The latest measurement data of
* either radar or laser.
*/
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Make sure you switch between lidar and radar
    measurements.
    */
    if (!is_initialized_) {

        // First measurement!

        cout << "EKF: " << endl;
        
        //Tune the inital values of state and covariance matricies:

        x_ << 1, 1, 1, 1, 1; // x & y position to be updated below

        P_ <<   std_laspx_, 0, 0, 0, 0, 
                0, std_laspy_, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, 1, 0,
                0, 0, 0, 0, 1;


        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            float rho = meas_package.raw_measurements_(0);
            float phi = meas_package.raw_measurements_(1);

            //Initialize State vector
            x_(0) = rho * cos(phi);
            x_(1) = rho * sin(phi);

            /*
            Note: Radar and CTRV velocity are not the same,
            so I cannot set ekf_.x_(2) = raw_measurements_(2); 
            */
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);
        }

        time_us_ = meas_package.timestamp_; // set the timestep for iteration 1

        //Done initializing, no need to predict or update
        is_initialized_ = true;
        return;

    }

    // Only Predict + Update for the sensor type if it is enabled
    if( (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) ||
        (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) ){

        /*
            Prediction Step
        */
        //compute the time elapsed between the current and previous measurements
        float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;    // expressed in seconds
        time_us_ = meas_package.timestamp_;

        //Call Predit Function
        Prediction(delta_t);

        /*
            Update Step
        */
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            //Radar
            UpdateRadar(meas_package);

        } else {
            //Laser
            UpdateLidar(meas_package);
        }

    }
}

/**
* Predicts sigma points, the state, and the state covariance matrix.
* @param {double} delta_t the change in time (in seconds) between the last
* measurement and this one.
*/
void UKF::Prediction(double delta_t) {
    /* 
        Generate Augmented Sigma Points
    */

    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1); // Init augmented sigma points
    VectorXd x_aug = VectorXd(n_aug_); // Init augmented state vector    
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_); //Init augmented covariance
        
    //set lambda for augmented sigma points
    lambda_ = 3 - n_aug_;

    // Create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    // Create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;

    // Create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    // Create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }    

    /* 
        Predict Augmented Sigma Points
    */
    
    for(int i = 0; i<2*n_aug_ +1; i++) // Iterate over augmented sigma points
    {
        //De-clutter:
        double p_x      = Xsig_aug(0, i);
        double p_y      = Xsig_aug(1, i);
        double v        = Xsig_aug(2, i);
        double yaw      = Xsig_aug(3, i);
        double yawd     = Xsig_aug(4, i);
        double nu_a     = Xsig_aug(5, i);
        double nu_yawdd = Xsig_aug(6, i);

        //Init Predictions:
        double px_p,py_p,v_p,yaw_p,yawd_p;

        //Calcuate the position state variables
        if(fabs(yawd) > 0.001) //Avoid divide by zero
        {
            //Compute positions:
            px_p = p_x + (v/yawd)*(sin(yaw + yawd*delta_t) - sin(yaw)) +
                    0.5*delta_t*delta_t*cos(yaw)*nu_a;
            py_p = p_y + (v/yawd)*(cos(yaw) - cos(yaw + yawd*delta_t)) +
                    0.5*delta_t*delta_t*sin(yaw)*nu_a;
        }
        else{
            // Assume moving straight
            px_p = p_x + v*delta_t * cos(yaw) +
                    0.5*delta_t*delta_t*cos(yaw)*nu_a;

            py_p = p_y + v*delta_t * sin(yaw) +
                    0.5*delta_t*delta_t*sin(yaw)*nu_a;
        }

        //Calculate other three state variables:
        v_p = v + 0 + delta_t * nu_a;
        yaw_p = yaw + yawd*delta_t + 0.5*delta_t*delta_t*nu_yawdd;
        yawd_p = yawd + 0 + delta_t * nu_yawdd;
        
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;

    }
    
    /* 
        Predict state and state covariance matrix
    */
    //Generate weights:
    weights_(0) = lambda_/(lambda_ + n_aug_);
    for(int i = 1; i<2*n_aug_ + 1; i++)
    {
        weights_(i) = 0.5/(n_aug_+lambda_);
    }

    //Predict State:
    x_.fill(0.0);
    for (int i = 0; i<2*n_aug_+1; i++) // Iterate over predicted sigma points
    {
        x_ = x_ + weights_(i)*Xsig_pred_.col(i);
    }

    //Predit Covariance: 
    P_.fill(0.0);   
    for (int i = 0; i<2*n_aug_+1; i++) // Iterate over predicted sigma points
    {        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        // angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();        
    }
}

/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */


    VectorXd z = meas_package.raw_measurements_; //Store the measurment vector

    int n_z = 2 ; // Define number of lidar points

    /* 
        Predict Measurment Step
    */

    // Create matrix for augmented sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    
    //Create matrix for predicted measurment model
    MatrixXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);

    // Translate predicted state vector into measurment space :
    for (int i = 0 ; i < 2 * n_aug_ + 1; i++)
    {
        Zsig(0,i) = Xsig_pred_(0,i);
        Zsig(1,i) = Xsig_pred_(1,i);
    }
 
    //Calculate Mean Predicted Measurement:
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //Calculate Measurement Covariance Matrix S:
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) 
    {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    } 

    // Add measurement noise to covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_laspx_*std_laspx_, 0, 
            0, std_laspy_*std_laspy_;
    S = S + R;


    /* 
        Perform UKF Update
    */


    //Calculate Cross Correlation Matrix 

    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) 
    {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    // Calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    // Update state mean and covariance matricies
    VectorXd z_diff = z - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    /* 
        Calculate NIS
    */

    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateRadar(MeasurementPackage meas_package) {

    VectorXd z = meas_package.raw_measurements_; //Store the measurment vector

    int n_z = 3 ; // Define number of radar points

    /* 
        Predict Measurment Step
    */

    // Create matrix for augmented sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    
    //Create matrix for predicted measurment model
    MatrixXd z_pred = VectorXd(n_z);
    
    // Translate predicted state vector into measurment space :
    for (int i = 0 ; i < 2 * n_aug_ + 1; i++)
    {
        // De-clutter
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;

        //Convert prediction vector to polar coordinates for radar measurment
        Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y); //r
        Zsig(1, i) = atan2(p_y, p_x); //phi
        Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
 
    //Calculate Mean Predicted Measurement:
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //Calculate Measurement Covariance Matrix S:
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) 
    {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    } 

    // Add measurement noise to covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0,std_radrd_*std_radrd_;
    S = S + R;

    /* 
        Perform UKF Update
    */

    //Calculate Cross Correlation Matrix T

    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) 
    {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    // Calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    // Update state mean and covariance matricies
    VectorXd z_diff = z - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    /* 
        Calculate NIS
    */

    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;



}

#include "ukf.h"
#include <Eigen/Dense>
#include <limits>

// debug
#include <iostream>
using std::cout;
using std::endl;

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Initializes Unscented Kalman Filter
 */
UKF::UKF()
{
    time_us_ = 0;
    is_initialized_ = false;


    // if this is false, laser measurement will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurement will be ignored (except during init)
    use_radar_ = true; 

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);


    /*
     * Measurement noise values
     * provided by the sensor manufacturer
     * NOT TO BE MODIFIED
     */
    
    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;
    
    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

    /// Variables to be initialized
    n_x_ = 5;
    n_aug_ = 7;
    n_sig_ = 2*n_aug_+1;
    lambda_ = 3 - n_aug_;

    // set weights
    weights_ = VectorXd(2*n_aug_ + 1);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    double weight = 0.5/(lambda_+n_aug_);
    for (int i=1; i<2*n_aug_+1; ++i)
        weights_(i) = weight;

    I_ = MatrixXd::Identity(5, 5);
    x_.fill(0.0);

    /// Variables to be tuned
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 3;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = M_PI/2.0;

    P_ = MatrixXd::Identity(n_x_, n_x_);
};

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
    if (!is_initialized_)
    {
        if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) 
        {
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);
        }
        is_initialized_ = true;
        time_us_ = meas_package.timestamp_;
        return;
    }
    
    double dt = (double)(meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;


    // Predict the state
    Prediction(dt);

    // Update the state
    if (use_laser_ && meas_package.sensor_type_ == meas_package.LASER)
    {
        // Lidar Measurement
        UpdateLidar(meas_package);
    }
    else if (use_radar_ && meas_package.sensor_type_ == meas_package.RADAR)
    {
        // Radar Measurement
        UpdateRadar(meas_package);
    }
}

void UKF::Prediction(double delta_t)
{
    // Create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);

    // Create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    // Create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

    // Create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0.0;
    x_aug(6) = 0.0;

    // create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5, 5) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;

    // Create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    // create augmented sigma points
    Xsig_aug.col(0) = x_aug;
    for (int i=0; i<n_aug_; ++i)
    {
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug.col(n_aug_ + i+1) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }

    double delta_t2 = delta_t * delta_t;

    // Predict sigma points
    Xsig_pred_ = MatrixXd(n_x_, Xsig_aug.cols());
    for (int i = 0; i < n_sig_; ++i)
    {
        double px = Xsig_aug(0, i);
        double py = Xsig_aug(1, i);
        double v = Xsig_aug(2, i);
        double yaw = Xsig_aug(3, i);
        double yawd = Xsig_aug(4, i);
        double nu_a = Xsig_aug(5, i);
        double nu_yawdd = Xsig_aug(6, i);

        // check division by zero
        if (fabs(yawd) > 0.001)
        {
            Xsig_pred_.col(i) << px + (v/yawd)*(sin(yaw+yawd*delta_t) - sin(yaw)) + (0.5*delta_t2*cos(yaw)*nu_a),
                                 py + (v/yawd)*(-cos(yaw+yawd*delta_t) + cos(yaw)) + (0.5*delta_t2*sin(yaw)*nu_a),
                                 v + (delta_t * nu_a),
                                 yaw + yawd*delta_t + (0.5*delta_t2*nu_yawdd),
                                 yawd + (delta_t*nu_yawdd);
        }
        else
        {
            Xsig_pred_.col(i) << px + v*cos(yaw)*delta_t + (0.5*delta_t2*cos(yaw)*nu_a),
                                 py + v*sin(yaw)*delta_t + (0.5*delta_t2*sin(yaw)*nu_a),
                                 v + (delta_t*nu_a),
                                 yaw + yawd*delta_t + (0.5*delta_t2*nu_yawdd),
                                 yawd + (delta_t*nu_yawdd);
        }
    }

    // predict state mean
    x_ = Xsig_pred_ * weights_;

    // Predict state covariance matrix
    P_.fill(0.0);
    for (int i=0; i<Xsig_pred_.cols(); ++i)
    {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        while (x_diff(3) >  M_PI) x_diff(3) -= 2.0 * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2.0 * M_PI;

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }
}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
    MatrixXd H(2, 5);
    H << 1, 0, 0, 0, 0,
         0, 1, 0, 0, 0;

    MatrixXd R(2, 2);
    R << std_laspx_*std_laspx_, 0.0,
         0.0, std_laspy_ * std_laspy_;
    
    // calculate innovation
    VectorXd y = meas_package.raw_measurements_ - (H*x_);
    MatrixXd Ht = H.transpose();
    MatrixXd S = H*P_*Ht + R;
    MatrixXd S_inv = S.inverse();
    
    MatrixXd K = (P_*Ht)*S_inv;

    x_ = x_ + K*y;
    P_ = (I_ - K*H)*P_;

    VectorXd nis = y.transpose() * S_inv * y;
    lidar_nis.push_back(nis(0));
}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
    int n_z = 3;
    // create matrix for sigma point in measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sig_);

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    
    // transform sigma points into measurement space 
    for (int i = 0; i < Xsig_pred_.cols(); ++i)
    {
        double px = Xsig_pred_(0, i);
        double py = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);

        Zsig(0, i) = sqrt(px*px + py*py);                             // radius
        Zsig(1, i)= atan2(py, px);                                    // phi angle
        Zsig(2, i) = (px*cos(yaw)*v + py*sin(yaw)*v)/Zsig(0, i); // radius_dot (velocity)
    }

    // Calculate mean predicted measurement
    z_pred.fill(0.0);
    z_pred = Zsig * weights_;

    // Calculate innovation covariance matrix S
    MatrixXd R(n_z, n_z); // Measurement sensor noise
    R.fill(0.0);
    R(0, 0) = std_radr_ * std_radr_;
    R(1, 1) = std_radphi_ * std_radphi_;
    R(2, 2) = std_radr_ * std_radr_;

    S.fill(0.0);
    for (int i = 0; i < n_sig_; ++i)
    {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    S = S + R;

    // Create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    // Calculate cross correlation Tc
    Tc.fill(0.0);
    for (int i = 0; i < n_sig_; ++i)
    {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

        VectorXd z_diff = Zsig.col(i) - z_pred;
        while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
       
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    // Calculate Kalman gain
    MatrixXd K = Tc * S.inverse();

    // Residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

    // angle normalization
    while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    x_ = x_ + K*z_diff;
    P_ = P_ - K*S*K.transpose();

    VectorXd nis = z_diff.transpose() * S.inverse() * z_diff;
    radar_nis.push_back(nis(0));
}

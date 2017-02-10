#pragma once

#include "Eigen/Eigen"
#include "Eigen/SVD"




class Alined{
public:

  enum DLT_METHOD{LINE_DLT,COMBINED_LINES} method_;

  Alined(DLT_METHOD = COMBINED_LINES);
  ~Alined();


  /*!
   * \brief Camera pose from line correspondences using DLT-Combined-Lines method
   *        Pose Estimation from Line Correspondences using Direct Linear Transformation" by
   *        Pribyl, B., Zemcik, P. and Cadik, M.
   * \param x - 3x(2N) 2D line endpoints (Don't need to correspond to 3D line endpoint locations)
   * \param X - 4x(2N) 3D line endpoints
   * \return
   */
  Eigen::Matrix4d poseFromLines(Eigen::Matrix<double,3,Eigen::Dynamic> x_c, Eigen::Matrix<double,4,Eigen::Dynamic> X_w);


private:

  /*!
   * \brief Create Plucker Line from two points only
   * \param X - 4x(2N) homogeneous point pair
   * \return L - 6x1 Plucker Line
   */
  Eigen::Matrix<double,6,Eigen::Dynamic> createPluckerLines(const Eigen::Matrix<double,4,Eigen::Dynamic> &X);

  /*!
   * \brief Calculate the Kronecker product with matrices m1 and m2
   * \param m1 - matrix
   * \param m2 - matrix
   * \return Matrix
   */
  Eigen::MatrixXd kron(Eigen::MatrixXd m1, Eigen::MatrixXd m2);

  /*!
   * \brief Create a skew-symmetric matrix from a vector
   * \param vec - 3x1 vector
   * \return skew 3x3 matrix
   */
  inline Eigen::Matrix3d skew(const Eigen::Vector3d& vec);

  /*!
   * \brief Create a vector from a skew-symmetric matrix
   * \param skew - 3x3 matrix
   * \return vec - 3x1 vector
   */
  inline Eigen::Vector3d unskew(const Eigen::Matrix3d& skew);

};



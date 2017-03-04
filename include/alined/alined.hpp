/*
* Copyright (C) 2016 Andrea Luca Lampart <lamparta at student dot ethz dot ch> (ETH Zurich)
* For more information see <https://github.com/andrealampart/alined>
*
* ALineD is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ALineD is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ALineD. If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include "Eigen/Eigen"
#include "Eigen/Geometry"
#include "Eigen/SVD"

#include <functional>


#define AL_LINE_DLT                 1
#define AL_COMBINED_LINES           2
#define AL_LEAST_SQUARES            4
#define AL_LEVENBERG_MARQUARDT      8
#define AL_NO_REFINE                16
#define AL_USE_REFINE               32
#define AL_HUBER_LOSS               64
#define AL_CAUCHY_LOSS              128

#define AL_MIN_LINES_COMBINED 5
#define AL_MIN_LINES_DLT 6
#define AL_MAX_ITER 1000

class Alined{
public:

  unsigned char method_;
  unsigned char solver_;
  unsigned char iterative_;

  Alined(unsigned char config);
  ~Alined();


  /*!
   * \brief Camera pose from line correspondences using DLT-Combined-Lines method
   *        Pose Estimation from Line Correspondences using Direct Linear Transformation" by
   *        Pribyl, B., Zemcik, P. and Cadik, M.
   * \param x_c - 3x(2N) 2D line endpoints (Don't need to correspond to 3D line endpoint locations)
   * \param X_w - 4x(2N) 3D line endpoints
   * \return Pose
   */
  Eigen::Matrix4d poseFromLines(Eigen::Matrix<double,3,Eigen::Dynamic> x_c, Eigen::Matrix<double,4,Eigen::Dynamic> X_w);

  /*!
   * \brief Camera pose using the iterative approach by Kumar and Hanson. This method needs a good prior.
   * \param x_c - 3x(2N) 2D line endpoints (Don't need to correspond to 3D line endpoint locations)
   * \param X_w - 4x(2N) 3D line endpoints
   * \return Pose
   */
  Eigen::Matrix4d poseFromLinesIterative(Eigen::Matrix4d pose, Eigen::Matrix<double,3,Eigen::Dynamic> x_c, Eigen::Matrix<double,4,Eigen::Dynamic> X_w);

  /*!
   * \brief Set scale of loss function
   * \param scale
   */
  void setLossScale(const double &scale);

private:

  double loss_scale_;

  /*!
   * \brief Wrapping any loss function
   */
  std::function<Eigen::Vector3d(const double&)> lossFunc_;


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

  /*!
   * \brief Iteratively find the correct pose using the R_and_T algorithm by Kumar and Hanson 1994
   * \param tf - Initial Pose
   * \param X - Point matrix
   * \param l_c - 2D line Matrix
   * \return Pose
   */
  Eigen::Matrix4d refineIteratively(const Eigen::Matrix4d &tf, Eigen::Matrix<double,4, Eigen::Dynamic> X_w, Eigen::Matrix<double, 3, Eigen::Dynamic> l_c, Eigen::Matrix<double, 1, Eigen::Dynamic> w);

  /*!
   * \brief Iteratively find the correct pose using the R_and_T algorithm by Kumar and Hanson 1994
   *        in a full Levenberg-Marquardt scheme. This reduces divergent behavior at close to singular situations.
   * \param tf - Initial Pose
   * \param X - Point matrix
   * \param l_c - 2D line Matrix
   * \return Pose
   */
  Eigen::Matrix4d levenbergMarquardt(const Eigen::Matrix4d &tf, Eigen::Matrix<double,4, Eigen::Dynamic> X_w, Eigen::Matrix<double, 3, Eigen::Dynamic> l_c, Eigen::Matrix<double, 1, Eigen::Dynamic> w);

  /*!
   * \brief Calculate the Huber loss function to penalize outliers in the nonlinear optimization
   * \param cost - incremental cost of measurement i
   * \param scale - cutoff-cost
   * \return weights - weigth to be used in outlier rejection
   */
  Eigen::Vector3d huberLoss(const double &cost);

  /*!
   * \brief Calculate the Cauchy loss function to penalize outliers in the nonlinear optimization
   * \param cost - incremental cost of measurement i
   * \param scale - cutoff-cost
   * \return weights - weigth to be used in outlier rejection
   */
  Eigen::Vector3d cauchyLoss(const double &cost);


};



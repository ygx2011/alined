#ifndef ALINED_HPP
#define ALINED_HPP

#include "Eigen/Eigen"




class Alined{

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


};

#endif // ALINED_HPP

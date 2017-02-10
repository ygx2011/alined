#include <ros/ros.h>
#include "alined/alined.hpp"
#include "geometry_msgs/PoseArray.h"


int main(int argc, char **argv)
{
  // Set up ROS.
  ros::init(argc, argv, "alined_test_node");

  std::cout << "Using Eigen v."<<EIGEN_WORLD_VERSION<<"."<<EIGEN_MAJOR_VERSION<<"."<<EIGEN_MINOR_VERSION<<"\n";

  Alined alined(alined.LINE_DLT);
  Eigen::MatrixXd x_c;
  Eigen::Matrix<double,4, 16> X_w;

  // Virtual camera position

  Eigen::Vector3d cam_pos_w(1.5, 0 , 0.5);
  Eigen::Matrix3d cam_ori_w;

  // Align camera's z-axis with world negative Y-axis
  // The choice of coordinate systems is consistent with the one in
  // "Pose Estimation from Line Correspondences using DLT" by B. Pribyl, 2016
  cam_ori_w << 1 ,0 ,0, 0, 0, -1, 0, 1, 0;

  Eigen::Matrix<double, 3, 4> projection_matrix;
  projection_matrix.block<3,3>(0,0) = cam_ori_w.inverse();
  projection_matrix.block<3,1>(0,3) = -(cam_ori_w.inverse())*cam_pos_w;

  //std::cout << "Projection Matrix = \n\n" << projection_matrix << "\n\n";

  //------------ Build 3D House ------------------//

  // (1,1,0)->(2,1,0)
  X_w.block<4,1>(0,0) = Eigen::Vector4d(0.2,1,0,1);
  X_w.block<4,1>(0,1) = Eigen::Vector4d(2,1,0,1);

  // (2,1,0)->(2,2,0)
  X_w.block<4,1>(0,2) = Eigen::Vector4d(1.1,1,1,1);
  X_w.block<4,1>(0,3) = Eigen::Vector4d(2.1,1,1,1);

  // (2,2,0)->(1,2,0)
  X_w.block<4,1>(0,4) = Eigen::Vector4d(1.2,1,2,1);
  X_w.block<4,1>(0,5) = Eigen::Vector4d(2.2,1,2,1);

  // (1,2,0)->(1,1,0)
  X_w.block<4,1>(0,6) = Eigen::Vector4d(1.12,1,1,1);
  X_w.block<4,1>(0,7) = Eigen::Vector4d(1.12,1,2,1);

  // (1,1,0)->(1,1,1)
  X_w.block<4,1>(0,8) = Eigen::Vector4d(0.98,1,0,1);
  X_w.block<4,1>(0,9) = Eigen::Vector4d(0.98,1,1,1);

  // (1,1,1)->(2,1,1)
  X_w.block<4,1>(0,10) = Eigen::Vector4d(6,1,2,1);
  X_w.block<4,1>(0,11) = Eigen::Vector4d(2,1,1,1);

  // (2,1,1)->(2,2,1)
  X_w.block<4,1>(0,12) = Eigen::Vector4d(1.9,1,1,1);
  X_w.block<4,1>(0,13) = Eigen::Vector4d(2,2,1,1);

  // (2,2,1)->(1,2,1)
  X_w.block<4,1>(0,14) = Eigen::Vector4d(2,2,1,1);
  X_w.block<4,1>(0,15) = Eigen::Vector4d(1,2,1,1);






  x_c = projection_matrix*X_w;


  std::cout <<"Begin Test: \n\n";
  alined.poseFromLines(x_c,X_w);

  std::cout << "Real Pose = \n\n"<< projection_matrix<<"\n\n";

  return 0;
}

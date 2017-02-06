#include <ros/ros.h>
#include "alined/alined.hpp"


int main(int argc, char **argv)
{
  // Set up ROS.
  ros::init(argc, argv, "ceres_test_node");

  Alined alined;
  Eigen::Matrix<double,3,10> x_c;
  Eigen::Matrix<double,4, 10> X_w;

  // Line 1
  X_w.block<4,1>(0,0) = Eigen::Vector4d(0,0,0,1);
  X_w.block<4,1>(0,1) = Eigen::Vector4d(1,0,0,1);

  x_c.block<3,1>(0,0) = Eigen::Vector3d(0,0,0);
  x_c.block<3,1>(1,1) = Eigen::Vector3d(0,0,0);

  X_w.block<4,1>(0,2) = Eigen::Vector4d(1,0,0,1);
  X_w.block<4,1>(0,3) = Eigen::Vector4d(2,0,0,1);

  x_c.block<3,1>(0,2) = Eigen::Vector3d(0,0,0);
  x_c.block<3,1>(1,3) = Eigen::Vector3d(0,0,0);

  X_w.block<4,1>(0,4) = Eigen::Vector4d(2,0,0,1);
  X_w.block<4,1>(0,5) = Eigen::Vector4d(3,0,0,1);

  x_c.block<3,1>(0,4) = Eigen::Vector3d(0,0,0);
  x_c.block<3,1>(1,5) = Eigen::Vector3d(0,0,0);

  X_w.block<4,1>(0,6) = Eigen::Vector4d(3,0,0,1);
  X_w.block<4,1>(0,7) = Eigen::Vector4d(4,0,0,1);

  x_c.block<3,1>(0,6) = Eigen::Vector3d(0,0,0);
  x_c.block<3,1>(1,7) = Eigen::Vector3d(0,0,0);

  X_w.block<4,1>(0,8) = Eigen::Vector4d(4,0,0,1);
  X_w.block<4,1>(0,9) = Eigen::Vector4d(5,0,0,1);

  x_c.block<3,1>(0,8) = Eigen::Vector3d(0,0,0);
  x_c.block<3,1>(1,9) = Eigen::Vector3d(0,0,0);


  //std::cout << X_w<< "\n\n" << X_w(4,0)<<"\n\n";

  std::cout <<"Begin Test: \n";
  alined.poseFromLines(x_c,X_w);

  return 0;
}

#include <ros/ros.h>
#include "alined/alined.hpp"
#include "geometry_msgs/PoseArray.h"
#include "external/tic_toc.hpp"


int main(int argc, char **argv)
{
  // Set up ROS.
  ros::init(argc, argv, "alined_test_node");


  Alined alined(alined.LINE_DLT, alined.USE_ITERATIVE_REFINEMENT);
  Eigen::MatrixXd x_c;
  Eigen::Matrix<double,4, 10> X_w;

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
  X_w.block<4,1>(0,0) = Eigen::Vector4d(0.02,1,0.0,1);
  X_w.block<4,1>(0,1) = Eigen::Vector4d(1.02,1,0.0,1);

  // (2,1,0)->(2,2,0)
  X_w.block<4,1>(0,2) = Eigen::Vector4d(0.0,1,0.21,1);
  X_w.block<4,1>(0,3) = Eigen::Vector4d(1.0,1,0.21,1);

  // (2,2,0)->(1,2,0)
  X_w.block<4,1>(0,4) = Eigen::Vector4d(0.0,1,0.43,1);
  X_w.block<4,1>(0,5) = Eigen::Vector4d(0.98,1,0.43,1);

  // (1,2,0)->(1,1,0)
  X_w.block<4,1>(0,6) = Eigen::Vector4d(0.04,1,0.62,1);
  X_w.block<4,1>(0,7) = Eigen::Vector4d(1.04,1,0.62,1);

  // (1,1,0)->(1,1,1)
  X_w.block<4,1>(0,8) = Eigen::Vector4d(0,1,0.8,1);
  X_w.block<4,1>(0,9) = Eigen::Vector4d(0.99,1,0.8,1);
/*
  // (1,1,1)->(2,1,1)
  X_w.block<4,1>(0,10) = Eigen::Vector4d(6,1,6,1);
  X_w.block<4,1>(0,11) = Eigen::Vector4d(2,1,1,1);

 // (2,1,1)->(2,2,1)
  X_w.block<4,1>(0,12) = Eigen::Vector4d(1.9,1,1,1);
  X_w.block<4,1>(0,13) = Eigen::Vector4d(2,2,1,1);

  // (2,2,1)->(1,2,1)
  X_w.block<4,1>(0,14) = Eigen::Vector4d(2,2,1,1);
  X_w.block<4,1>(0,15) = Eigen::Vector4d(1,2,1,1);*/






  x_c = projection_matrix*X_w;
   Eigen::Matrix<double, 3, Eigen::Dynamic> noise;
   noise.setRandom(3,x_c.cols());
   noise = noise.eval()*0.01-Eigen::MatrixXd().setOnes(3,x_c.cols())*0.005;

   std::cout << noise<<"\n\n";


  std::cout <<"Test 1: \n\n";
  //Eigen::Matrix4d tf = alined.poseFromLines(x_c,X_w);



  std::cout << "Real Pose = \n\n"<< projection_matrix<<"\n\n";
  //std::cout << "Estimated Pose = \n\n"<< tf <<"\n\n";

  std::cout <<"Test 2: \n\n";


  Eigen::Matrix4d tf_2;
  Eigen::Matrix3d disturbance(Eigen::AngleAxisd(0, Eigen::Vector3d(0,0,1)));

  tf_2.block<3,3>(0,0) = disturbance*projection_matrix.block<3,3>(0,0);
  tf_2.block<3,1>(0,3) = projection_matrix.block<3,1>(0,3);
  tf_2.block<1,4>(3,0) = Eigen::Vector4d(0,0,0,1).transpose();

  for(int i = 0; i<10; i++){

    projection_matrix.block<3,1>(0,3) = projection_matrix.block<3,1>(0,3).eval() -(cam_ori_w.inverse())*Eigen::Vector3d(0,-1,0.01);
    x_c = projection_matrix*X_w;
    x_c = x_c.eval() + noise;
    tic("nsec","Iterative");
    tf_2 = alined.poseFromLinesIterative(tf_2,x_c,X_w);
    toc();
    std::cout << "tf(" << i << ") = \n\n"<< tf_2 <<"\n\n";

  }

  double err_rot = acos((tf_2.block<3,3>(0,0).transpose()*projection_matrix.block<3,3>(0,0)).trace()/2-0.5);
  double err_trans = (tf_2.block<3,1>(0,3)-projection_matrix.block<3,1>(0,3)).norm();
  double noise_level = (noise.array()/x_c.array()).sum()/noise.cols();
  std::cout << "Lines = " << x_c.cols()/2 << "  Rotational Error = " << err_rot*180/3.14158 << "  Translation error = "<<  err_trans << "  Relative Noise = "<<noise_level<<"\n";

  return 0;
}

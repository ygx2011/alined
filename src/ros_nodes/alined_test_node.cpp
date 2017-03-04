#include <ros/ros.h>
#include "alined/alined.hpp"
#include "geometry_msgs/PoseStamped.h"
#include "external/tic_toc.hpp"
#include "visualization_msgs/MarkerArray.h"


int main(int argc, char **argv)
{
  // Set up ROS.
  ros::init(argc, argv, "alined_test_node");
  ros::NodeHandle nh;

  const int nlines = 4;

  Alined alined(AL_COMBINED_LINES|AL_NO_REFINE|AL_LEVENBERG_MARQUARDT|AL_CAUCHY_LOSS);
  alined.setLossScale(1.0);

  Eigen::MatrixXd x_c;
  Eigen::Matrix<double,4, 2*nlines> X_w;
  X_w.setRandom(4,2*nlines);
  X_w.block<1,2*nlines>(3,0) = Eigen::MatrixXd::Ones(1,2*nlines);

  // Virtual camera position

  Eigen::Vector3d cam_pos_w(0, 0 , 2);
  Eigen::Matrix3d cam_ori_w;

  // Align camera's z-axis with world negative Y-axis
  // The choice of coordinate systems is consistent with the one in
  // "Pose Estimation from Line Correspondences using DLT" by B. Pribyl, 2016
  cam_ori_w << 1 ,0 ,0, 0, 0, -1, 0, 1, 0;

  Eigen::Matrix<double, 3, 4> projection_matrix;
  projection_matrix.block<3,3>(0,0) = cam_ori_w.inverse();
  projection_matrix.block<3,1>(0,3) = -(cam_ori_w.inverse())*cam_pos_w;

  //std::cout << "Projection Matrix = \n\n" << projection_matrix << "\n\n";

  //------------ Build 3D shelf scene ------------------//

  // (1,1,0)->(2,1,0)
  /*X_w.block<4,1>(0,0) = Eigen::Vector4d(0.02,1,0.0,1);
  X_w.block<4,1>(0,1) = Eigen::Vector4d(1.02,1,0.0,1);

  // (2,1,0)->(2,2,0)
  X_w.block<4,1>(0,2) = Eigen::Vector4d(0.02,1,0.0,1);
  X_w.block<4,1>(0,3) = Eigen::Vector4d(1.02,1,0.0,1);

  // (2,2,0)->(1,2,0)
  X_w.block<4,1>(0,4) = Eigen::Vector4d(0.0,1,0.43,1);
  X_w.block<4,1>(0,5) = Eigen::Vector4d(0.98,1,0.43,1);

  // (1,2,0)->(1,1,0)
  X_w.block<4,1>(0,6) = Eigen::Vector4d(0.04,1,0.62,1);
  X_w.block<4,1>(0,7) = Eigen::Vector4d(1.04,1,0.62,1);

  // (1,1,0)->(1,1,1)
  X_w.block<4,1>(0,8) = Eigen::Vector4d(0,1,0.8,1);
  X_w.block<4,1>(0,9) = Eigen::Vector4d(0.99,1,0.8,1);

  // (1,1,1)->(2,1,1)
  X_w.block<4,1>(0,10) = Eigen::Vector4d(2.1,1,1.1,1);
  X_w.block<4,1>(0,11) = Eigen::Vector4d(2,1,1,1);

 // (2,1,1)->(2,2,1)
  X_w.block<4,1>(0,12) = Eigen::Vector4d(1.9,1,1,1);
  X_w.block<4,1>(0,13) = Eigen::Vector4d(2,1,1,1);

  // (2,2,1)->(1,2,1)
  X_w.block<4,1>(0,14) = Eigen::Vector4d(2,2,1,1);
  X_w.block<4,1>(0,15) = Eigen::Vector4d(1,2,1,1);
*/





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

  tf_2.block<3,3>(0,0) = projection_matrix.block<3,3>(0,0);
  tf_2.block<3,1>(0,3) = projection_matrix.block<3,1>(0,3);
  tf_2.block<1,4>(3,0) = Eigen::Vector4d(0,0,0,1).transpose();


  ros::Publisher posePub = nh.advertise<geometry_msgs::PoseStamped>("Pose", 10);
  ros::Publisher posePub_2 = nh.advertise<geometry_msgs::PoseStamped>("Pose_Orig", 10);

  ros::Publisher line_pub_ = nh.advertise<visualization_msgs::Marker>("lines",10);
  visualization_msgs::Marker lines;
  lines.header.frame_id = "world";
  lines.action = visualization_msgs::Marker::ADD;
  lines.type = visualization_msgs::Marker::LINE_LIST;
  std_msgs::ColorRGBA color;
  color.r = 0.0f;
  color.g = 1.0f;
  color.b = 0.0f;
  color.a = 1.0f;
  lines.color = color;

  lines.scale.x=0.005;




  for(int i=0; i< X_w.cols();++i){
    geometry_msgs::Point p;
    p.x = X_w(0,i);
    p.y = X_w(1,i);
    p.z = X_w(2,i);
    lines.points.push_back(p);
    lines.colors.push_back(color);
  }



  ros::Time::init();

  for(int i = 0; i<100; i++){
    line_pub_.publish(lines);

    Eigen::AngleAxisd a(0.3,Eigen::Vector3d(0,0,1));
    Eigen::Quaterniond dr(a);

    projection_matrix.block<3,3>(0,0) = dr.toRotationMatrix()*projection_matrix.block<3,3>(0,0).eval();

    projection_matrix.block<3,1>(0,3) = projection_matrix.block<3,1>(0,3).eval() -(projection_matrix.block<3,3>(0,0).inverse())*Eigen::Vector3d(0,-0.1,0.1);
    x_c = projection_matrix*X_w;
    x_c = x_c.eval()+ noise;
    tic("nsec","Iterative");
    tf_2 = alined.poseFromLinesIterative(tf_2,x_c,X_w);
    //tf_2 = alined.poseFromLines(x_c,X_w);
    toc();
    std::cout << "tf(" << i << ") = \n\n"<< tf_2 <<"\n\n" << projection_matrix << "\n\n";


    geometry_msgs::PoseStamped pose_stamped;
    pose_stamped.header.frame_id = "world";
    pose_stamped.header.stamp = ros::Time::now();

    geometry_msgs::Point point;
    point.x = tf_2(0,3);
    point.y = tf_2(1,3);
    point.z = tf_2(2,3);

    geometry_msgs::Quaternion quat;

    Eigen::Quaterniond r(tf_2.block<3,3>(0,0));

    quat.x = r.x();
    quat.y = r.y();
    quat.z = r.z();
    quat.w = r.w();

    pose_stamped.pose.position = point;

    pose_stamped.pose.orientation = quat;

    posePub.publish(pose_stamped);


    geometry_msgs::PoseStamped pose_stamped_2;
    pose_stamped_2.header.frame_id = "world";
    pose_stamped_2.header.stamp = ros::Time::now();


    point.x = projection_matrix(0,3);
    point.y = projection_matrix(1,3);
    point.z = projection_matrix(2,3);



    Eigen::Quaterniond r_2(projection_matrix.block<3,3>(0,0));

    quat.x = r_2.x();
    quat.y = r_2.y();
    quat.z = r_2.z();
    quat.w = r_2.w();

    pose_stamped_2.pose.position = point;

    pose_stamped_2.pose.orientation = quat;

    posePub_2.publish(pose_stamped_2);
    ros::Duration(0.05).sleep();

  }

  double err_rot = acos((tf_2.block<3,3>(0,0).transpose()*projection_matrix.block<3,3>(0,0)).trace()/2-0.5);
  double err_trans = (tf_2.block<3,1>(0,3)-projection_matrix.block<3,1>(0,3)).norm();
  double noise_level = (noise.array()/x_c.array()).sum()/noise.cols();
  std::cout << "Lines = " << x_c.cols()/2 << "  Rotational Error = " << err_rot*180/3.14158 << "  Translation error = "<<  err_trans << "  Relative Noise = "<<noise_level<<"\n";

  return 0;
}

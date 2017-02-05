#include "alined/alined.hpp"
#include "exception"
#include "iostream"

#define MIN_LINES 5

Eigen::Matrix4d Alined::poseFromLines(Eigen::Matrix<double,3,Eigen::Dynamic> x_c, Eigen::Matrix<double,4, Eigen::Dynamic> X_w){

// Following the matlab implementation of the original authors...

  // Check input size
  if(X_w.cols()%2 != 0){
    std::cout << "3D line endpoint vector needs to have even column size. Check dimesions...\n";
    throw std::exception();
  }

  if(x_c.cols()!=X_w.cols()){
    std::cout << "3D line endpoints and 2D line endpoints need to be of equal column size. Check dimensions...\n";
    throw std::exception();
  }

  int nlines = (int)(X_w.cols()/2);

  if(nlines < MIN_LINES){
    std::cout << "Not enough lines to extract pose from. Required: "<< MIN_LINES <<" Provided: "<< nlines <<"\n";
    throw std::exception();
  }

  // Create Plucker Lines from endpoints
  Eigen::Matrix<double,6,Eigen::Dynamic> L_w = createPluckerLines(X_w);
  
  // Prenormalization

  // Create line measurements vectors
  Eigen::Matrix<double, 3, Eigen::Dynamic> x_1 = x_c.Map(x_c.data(),6, x_c.cols()/2).topRows(3);
  Eigen::Matrix<double, 3, Eigen::Dynamic> x_2 = x_c.Map(x_c.data(),6, x_c.cols()/2).bottomRows(3);

  //todo: cross product of each
  //l_c =

  // Combined Measurement Matrix




}


Eigen::Matrix<double,6,Eigen::Dynamic> createPluckerLines(const Eigen::Matrix<double,4,Eigen::Dynamic> &X){

  Eigen::Matrix<double, 4, Eigen::Dynamic> X_1 = X.Map(X.data(),8, X.cols()/2).topRows(4);
  Eigen::Matrix<double, 4, Eigen::Dynamic> X_2 = X.Map(X.data(),8, X.cols()/2).bottomRows(4);


  Eigen::Matrix<double,6, Eigen::Dynamic> L;
  return L;

}

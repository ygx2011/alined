#include "alined/alined.hpp"
#include "exception"
#include "iostream"

#define MIN_LINES 5

Alined::Alined(){

}

Alined::~Alined(){

}

Eigen::Matrix4d Alined::poseFromLines(Eigen::Matrix<double,3,Eigen::Dynamic> x_c, Eigen::Matrix<double,4, Eigen::Dynamic> X_w){

// Following the matlab implementation of the original authors...

  //----------------------- Check input size ----------------------------------------------------------------------//
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

  //----------------- Create Plucker Lines from endpoints----------------------------------------------------------//

  Eigen::Matrix<double,6,Eigen::Dynamic> L_w = createPluckerLines(X_w);
  
  //----------------------- Prenormalization ----------------------------------------------------------------------//

  // Normalize the V part of the pluecker lines s.t. ||V|| = sqrt(3) [Bartoli, Sturm, 2004]

  for(int i = 0; i < nlines; i++){
    Eigen::Vector3d V = L_w.block<3,1>(3,i);
    L_w.block<6,1>(0,i) = (sqrt(3) * L_w.block<6,1>(0,i).eval())/V.norm();
  }




  //----------------- Create line measurements vectors-------------------------------------------------------------//

  Eigen::Matrix<double, 3, Eigen::Dynamic> x_1 = Eigen::MatrixXd::Map(x_c.data(),6, x_c.cols()/2).topRows(3);
  Eigen::Matrix<double, 3, Eigen::Dynamic> x_2 = Eigen::MatrixXd::Map(x_c.data(),6, x_c.cols()/2).bottomRows(3);

  Eigen::Matrix<double, 3, Eigen::Dynamic> l_c;
  l_c.resize(3, x_1.cols());

  for(int i = 0; i < nlines; i++){
    l_c.block<3,1>(0,i) = x_1.block<3,1>(0,i).eval().cross(x_2.block<3,1>(0,i).eval());
  }

  std::cout << "Measurement Matrix l_c = \n\n"<<l_c<<"\n\n";

  //---------------- Combined Measurement Matrix ------------------------------------------------------------------//

  Eigen::Matrix<double, 4, Eigen::Dynamic> X_1_w = Eigen::MatrixXd::Map(X_w.data(), 8, X_w.cols()/2).topRows(4);
  Eigen::Matrix<double, 4, Eigen::Dynamic> X_2_w = Eigen::MatrixXd::Map(X_w.data(), 8, X_w.cols()/2).bottomRows(4);

  // Measurement matrix for point - line correspondence



  /*Eigen::Matrix2d a;
  a << 1,2,3,4;
  Eigen::Matrix2d b;
  b << 1,1,1,1;

  std::cout<< a << "\n\n" << b << "\n \n" << kron(a,b) <<"\n\n";*/



}


Eigen::MatrixXd Alined::kron(Eigen::MatrixXd m1, Eigen::MatrixXd m2){

  Eigen::MatrixXd result(m1.rows()*m2.rows(), m1.cols()*m2.cols());

  for(int i = 0; i < m1.cols(); i++){
    for(int j = 0; j < m1.rows(); j++){
      result.block(i*m2.rows(),j*m2.cols(),m2.rows(),m2.cols()) = m1(j,i)*m2;
    }
  }

  return result;

}


Eigen::Matrix<double,6,Eigen::Dynamic> Alined::createPluckerLines(const Eigen::Matrix<double,4,Eigen::Dynamic> &X){

  /*
   * Decompose X into X_1 and X_2 through remapping
   */

  Eigen::Matrix<double, 4, Eigen::Dynamic> X_1 = Eigen::MatrixXd::Map(X.data(), 8,X.cols()/2).topRows(4);
  Eigen::Matrix<double, 4, Eigen::Dynamic> X_2 = Eigen::MatrixXd::Map(X.data(), 8, X.cols()/2).bottomRows(4);

  //std::cout << X_1 <<"\n\n" << X_2 << "\n\n";

  /* Plücker Line is created from Plücker matrix, such that (row-major)
   *
   *  L = {l12, l20, l01, l30,l31,l32} represents {X1 x X2; inhom(X2-X1)}
   *
   */

  // Need to compute L at runtime (Dynamic needed)
  Eigen::Matrix<double, 6, Eigen::Dynamic> L;
  L.resize(6,X_1.cols());

  for(int i = 0; i < X_1.cols(); i++){
    L.block<3,1>(0,i) = X_1.block<3,1>(0,i).eval().cross(X_2.block<3,1>(0,i).eval());
    L.block<3,1>(3,i) = X_2.block<3,1>(0,i).eval()-X_1.block<3,1>(0,i).eval();
  }

  std::cout << "3D Line Matrix = \n\n" << L << "\n\n";

  return L;

}

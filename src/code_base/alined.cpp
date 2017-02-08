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
  /* The practicability of the Prenormalization is explained in more detail in
   * "In Defense of the Eight-Point Algorithm" [Hartley 1997]
   *
   * Both translation and anisotropic scaling will be used.
   */

  // We assume each Point to have unit 4th component. No normalization is needed here.

  // Normalize the V part of the pluecker lines s.t. ||V|| = sqrt(3) [Bartoli, Sturm, 2004]
  // As the moment V is coupled to the distance from the line to the origin, this places the line
  // in [1 1 1] distance from the origin.

  /*for(int i = 0; i < nlines; i++){
    Eigen::Vector3d V = L_w.block<3,1>(3,i);
    L_w.block<6,1>(0,i) = (sqrt(3) * L_w.block<6,1>(0,i).eval())/V.norm();
  }*/

  // Translate points and lines to have their centroid at the origin





  //----------------- Create line measurements vectors-------------------------------------------------------------//

  Eigen::Matrix<double, 3, Eigen::Dynamic> x_1 = Eigen::MatrixXd::Map(x_c.data(),6, x_c.cols()/2).topRows(3);
  Eigen::Matrix<double, 3, Eigen::Dynamic> x_2 = Eigen::MatrixXd::Map(x_c.data(),6, x_c.cols()/2).bottomRows(3);

  Eigen::Matrix<double, 3, Eigen::Dynamic> l_c;
  l_c.resize(3, nlines);

  for(int i = 0; i < nlines; i++){
    l_c.block<3,1>(0,i) = x_1.block<3,1>(0,i).eval().cross(x_2.block<3,1>(0,i).eval());
  }

  //std::cout << "Measurement Matrix l_c = \n\n"<<l_c<<"\n\n";

  //---------------- Combined Measurement Matrix ------------------------------------------------------------------//

  Eigen::Matrix<double, 4, Eigen::Dynamic> X_1_w = Eigen::MatrixXd::Map(X_w.data(), 8, nlines).topRows(4);
  Eigen::Matrix<double, 4, Eigen::Dynamic> X_2_w = Eigen::MatrixXd::Map(X_w.data(), 8, nlines).bottomRows(4);


  // Measurement matrix for point - line correspondence
   Eigen::Matrix<double, 3, Eigen::Dynamic> l_c_l_c;
   l_c_l_c.resize( 3, 2 * nlines);
   l_c_l_c.block( 0, 0, 3,nlines) = l_c;
   l_c_l_c.block(0,nlines, 3, nlines) = l_c;

   Eigen::Matrix<double, 4, Eigen::Dynamic> X_1_X_2;
   X_1_X_2.resize(4,2*nlines);
   X_1_X_2.block(0,0,4,nlines) = X_1_w;
   X_1_X_2.block(0,nlines,4,nlines) = X_2_w;



  Eigen::MatrixXd M_pl_one = kron(Eigen::Vector4d(1.0, 1.0, 1.0, 1.0).transpose(), l_c_l_c.transpose());

  //std::cout << M_tmp_one <<"\n\n";

  Eigen::MatrixXd M_pl_two = kron(X_1_X_2.transpose(), Eigen::Vector3d(1.0,1.0,1.0).transpose());

  //std::cout << M_tmp_two <<"\n\n";

  Eigen::MatrixXd M_pl = (M_pl_one.array()*M_pl_two.array()).matrix();


  // Measurement matrix for line - line correspondence
  Eigen::MatrixXd six_ones;
  six_ones.setOnes(1,6);

  Eigen::MatrixXd tmp_mat;
  tmp_mat.setZero(nlines*2, 3);

  tmp_mat.block(0,0,nlines,1) = l_c.block(2,0,1,nlines).transpose();
  tmp_mat.block(0,2,nlines,1) = (-1)*l_c.block(0,0,1,nlines).transpose();
  tmp_mat.block(nlines,1,nlines,1) = l_c.block(2,0,1,nlines).transpose();
  tmp_mat.block(nlines,2,nlines,1) = (-1)*l_c.block(1,0,1,nlines).transpose();

  //std::cout <<"tmp = \n\n" << tmp_mat << "\n\n";

  Eigen::Vector3d three_ones(1,1,1);

  Eigen::MatrixXd l_w_l_w;
  l_w_l_w.resize(6,nlines*2);

  l_w_l_w.block(0,0,6,nlines) = L_w;
  l_w_l_w.block(0,nlines,6,nlines) = L_w;

  Eigen::MatrixXd M_ll_one = kron(six_ones,tmp_mat);
  Eigen::MatrixXd M_ll_two = kron(l_w_l_w.transpose(), three_ones.transpose());

  Eigen::MatrixXd M_ll = (M_ll_one.array()*M_ll_two.array()).matrix();

  // Normalize measurement matrices

  double M_pl_ss = (M_pl.array()*M_pl.array()).sum();
  double M_ll_ss = (M_ll.array()*M_ll.array()).sum();

  M_pl = M_pl.eval()/sqrt(M_pl_ss);
  M_ll = M_ll.eval()/sqrt(M_ll_ss);

  //std::cout << "M_ll = \n\n" << M_ll << "\n\n";
  //std::cout << "M_pl = \n\n" << M_pl << "\n\n";

  // Combine measurement matrices
  Eigen::MatrixXd M;
  M.setZero(M_pl.rows()+M_ll.rows(),21);

  M.block(0,0,M_pl.rows(),M_pl.cols()) = M_pl;
  M.block(M_pl.rows(),0, M_ll.rows(),9) = M_ll.block(0,0,M_ll.rows(),9);
  M.block(M_pl.rows(),12, M_ll.rows(),9) = M_ll.block(0,9,M_ll.rows(),9);

  //std::cout << M << "\n\n";

  //------------------------------------ Pose Estimation ------------------------//

  // Perform a singular value decomposition
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullV);
  Eigen::MatrixXd V_right = svd.matrixV().rightCols(1);
  V_right.resize(3,7);

  // First estimate of projection matrix using last singular vector
  Eigen::Matrix<double, 3, 7> P_est = V_right;

  // Divide combined projection matrix into submatrices
  Eigen::MatrixXd P1_est = P_est.leftCols(3);
  Eigen::MatrixXd P2_est = P_est.block<3,1>(0,3);
  Eigen::MatrixXd P3_est = P_est.rightCols(3);

  //std::cout << "P1_est = " << P1_est << "\n\n";

  // Algorithm 1 + 2 from paper: scale estimation + R1 orthogonalization
  Eigen::JacobiSVD<Eigen::MatrixXd> svd_P(P1_est, Eigen::ComputeFullV|Eigen::ComputeFullU);
  double det = (svd_P.matrixU()*((svd_P.matrixV()).transpose())).determinant();
  double scale = det/svd_P.singularValues().array().sum()*3;

  //std::cout << "U = " << svd_P.matrixU() <<"\n\n";
  //std::cout << "V = " << svd_P.matrixV() <<"\n\n";

  std::cout << "Condition Factor = " << ((svd_P.singularValues())(0,0))/(svd_P.singularValues()(2,0)) <<"\n\n";


  Eigen::MatrixXd R_est = det*(svd_P.matrixU()*svd_P.matrixV().transpose());
  Eigen::MatrixXd T_est = scale*R_est.transpose()*P2_est;


  std::cout << "R_est = \n\n" << R_est <<"\n\n" << "T_est = \n\n" << -T_est << "\n\n";




}


Eigen::MatrixXd Alined::kron(Eigen::MatrixXd m1, Eigen::MatrixXd m2){

  Eigen::MatrixXd result(m1.rows()*m2.rows(), m1.cols()*m2.cols());

  for(int i = 0; i < m1.cols(); i++){
    for(int j = 0; j < m1.rows(); j++){
      result.block(j*m2.rows(),i*m2.cols(),m2.rows(),m2.cols()) = m1(j,i)*m2;
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

  //std::cout << "3D Line Matrix = \n\n" << L << "\n\n";

  return L;

}

#include "alined/alined.hpp"
#include "exception"
#include "iostream"
#include "unsupported/Eigen/MatrixFunctions"
#include "external/tic_toc.hpp"

#define MIN_LINES_COMBINED 5
#define MIN_LINES_DLT 6

Alined::Alined(DLT_METHOD method)
  :method_(method)
{


}

Alined::~Alined(){

}

Eigen::Matrix4d Alined::poseFromLines(Eigen::Matrix<double,3,Eigen::Dynamic> x_c, Eigen::Matrix<double,4, Eigen::Dynamic> X_w){


  tic("nsec", "Pose Estimation:");

  //Output matrix
  Eigen::Matrix4d tf;

//----------------------------DLT-COMBINED-LINES-------------------------------------------------------------------//
  if(method_ == COMBINED_LINES){

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

  if(nlines < MIN_LINES_COMBINED){
    std::cout << "Not enough lines to extract pose from. Required: "<< MIN_LINES_COMBINED <<" Provided: "<< nlines <<"\n";
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

  for(int i = 0; i < nlines; i++){
    Eigen::Vector3d V = L_w.block<3,1>(3,i);
    L_w.block<6,1>(0,i) = (sqrt(3) * L_w.block<6,1>(0,i).eval())/V.norm();
  }

  // Translate points and lines to have their centroid at the origin
  double T_prenorm_x = X_w.block(0,0,1,X_w.cols()).array().sum()/X_w.cols();
  double T_prenorm_y = X_w.block(1,0,1,X_w.cols()).array().sum()/X_w.cols();
  double T_prenorm_z = X_w.block(2,0,1,X_w.cols()).array().sum()/X_w.cols();

  Eigen::Vector3d T_prenorm(T_prenorm_x,T_prenorm_y,T_prenorm_z);
  Eigen::Matrix4d D1_pl;
  D1_pl.block<3,3>(0,0) = Eigen::DiagonalMatrix<double,3,3>(1,1,1);
  D1_pl.block<3,1>(0,3) = T_prenorm;
  D1_pl.block<1,4>(3,0) = Eigen::Vector4d(0,0,0,1).transpose();

  Eigen::Matrix<double, 6,6> D1_ll;
  D1_ll.block<3,3>(0,0) = Eigen::DiagonalMatrix<double,3,3>(1,1,1);
  D1_ll.block<3,3>(0,3) = skew(T_prenorm);
  D1_ll.block<3,3>(3,0) = Eigen::Matrix3d::Zero();
  D1_ll.block<3,3>(3,3) = Eigen::DiagonalMatrix<double,3,3>(1,1,1);

  std::cout << D1_pl << "\n\n" << D1_ll << "\n\n";



  //----------------- Create line measurements vectors-------------------------------------------------------------//

  Eigen::Matrix<double, 3, Eigen::Dynamic> x_1 = Eigen::MatrixXd::Map(x_c.data(),6, x_c.cols()/2).topRows(3);
  Eigen::Matrix<double, 3, Eigen::Dynamic> x_2 = Eigen::MatrixXd::Map(x_c.data(),6, x_c.cols()/2).bottomRows(3);

  Eigen::Matrix<double, 3, Eigen::Dynamic> l_c;
  l_c.resize(3, nlines);

  for(int i = 0; i < nlines; i++){
    l_c.block<3,1>(0,i) = x_1.block<3,1>(0,i).cross(x_2.block<3,1>(0,i));
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
  double scale = 3*det/(svd_P.singularValues().array().sum());

  //std::cout << "U = " << svd_P.matrixU() <<"\n\n";
  //std::cout << "V = " << svd_P.matrixV() <<"\n\n";

  std::cout << "Condition Factor = " << ((svd_P.singularValues())(0,0))/(svd_P.singularValues()(2,0)) <<"\n\n";


  Eigen::MatrixXd R_est = det*(svd_P.matrixU()*svd_P.matrixV().transpose());
  Eigen::MatrixXd T_est = scale*R_est.transpose()*P2_est;


  std::cout << "R_est = \n\n" << R_est <<"\n\n" << "T_est = \n\n" << -T_est << "\n\n";

  // Algorithm 3: decomposition of an essential matrix R[t]x
  // from "Uniqueness and estimation of 3D motion...," , Tsai, Huang 1984

  Eigen::JacobiSVD<Eigen::MatrixXd> svd_P_3(scale*P3_est, Eigen::ComputeFullV|Eigen::ComputeFullU);

  Eigen::Matrix3d Z;
  Z << 0,1,0,-1,0,0,0,0,0;

  Eigen::Matrix3d W;
  W << 0,-1,0,1,0,0,0,0,1;
  double q = (svd_P_3.singularValues().block<2,1>(0,0).array().sum())*0.5;


  // 2 possible Solutions A/B
  double det_A = (svd_P_3.matrixU()*W*svd_P_3.matrixV().transpose()).determinant();
  double det_B = (svd_P_3.matrixU()*W.transpose()*svd_P_3.matrixV().transpose()).determinant();



  Eigen::Matrix3d R_A = svd_P_3.matrixU()*W*Eigen::DiagonalMatrix<double,3>(1,1,det_A)*svd_P_3.matrixV().transpose();
  Eigen::Matrix3d R_B = svd_P_3.matrixU()*W.transpose()*Eigen::DiagonalMatrix<double,3>(1,1,det_B)*svd_P_3.matrixV().transpose();

  Eigen::Matrix3d T_A = q*svd_P_3.matrixV()*Z*svd_P_3.matrixV().transpose();
  Eigen::Matrix3d T_B = q*svd_P_3.matrixV()*Z.transpose()*svd_P_3.matrixV().transpose();

  //Calculate nearest skew symmetric matrix in case of noisy values
  T_A = (T_A.eval()-T_A.eval().transpose())*0.5;
  T_B = (T_B.eval()-T_B.eval().transpose())*0.5;

  //std::cout << det_A << "\n\n" << det_B << "\n\n";


  std::cout << "R_A =\n\n" <<R_A <<"\n\nR_B =\n\n"<<R_B<<"\n\n";

  // Choose a solution based on which one is closer to R_est (closer to viewing direction)
  double theta_A = acos(((R_A.transpose()*R_est).trace()-1)*0.5);
  double theta_B = acos(((R_B.transpose()*R_est).trace()-1)*0.5);

  Eigen::Matrix3d R_est_two;
  Eigen::Vector3d T_est_two;

  if(theta_A <= theta_B){
    R_est_two = R_A;
    T_est_two = unskew(T_A);
  }
  else{
    R_est_two = R_B;
    T_est_two = unskew(T_B);
  }



  std::cout << "R_est_two =\n\n" <<R_est_two <<"\n\nT_est_two =\n\n"<<T_est_two<<"\n\n";

  // Generate output
  double k = 0.7;

  //Interpolate rotations

  Eigen::Matrix3d log_R_diff = (R_est_two.transpose()*R_est).log();
  Eigen::Matrix3d R = R_est*(((1-k)*log_R_diff).exp());
  Eigen::Vector3d T = ((k*T_est).array() + ((1-k)*T_est_two).array()).matrix();

  tf.block<3,3>(0,0) = R;
  tf.block<3,1>(0,3) = R*T;
  tf.block<1,4>(3,0) = Eigen::Vector4d(0,0,0,1).transpose();

  std::cout << "Estimated Pose = \n\n"<<tf<<"\n\n";

  }

  //----------------------------------DLT-LINES----------------------------------------------------------------------//
  else if(method_ == LINE_DLT){

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

    if(nlines < MIN_LINES_DLT){
      std::cout << "Not enough lines to extract pose from. Required: "<< MIN_LINES_DLT <<" Provided: "<< nlines <<"\n";
      throw std::exception();
    }

    //----------------------- Prenormalization ----------------------------------------------------------------------//
    /* The practicability of the Prenormalization is explained in more detail in
     * "In Defense of the Eight-Point Algorithm" [Hartley 1997]
     *
     * Both translation and anisotropic scaling will be used.
     */

    // We assume each Point to have unit 4th component. No normalization is needed here.

    // Translate points and lines to have their centroid at the origin
    double T_prenorm_x = X_w.block(0,0,1,X_w.cols()).array().sum()/X_w.cols();
    double T_prenorm_y = X_w.block(1,0,1,X_w.cols()).array().sum()/X_w.cols();
    double T_prenorm_z = X_w.block(2,0,1,X_w.cols()).array().sum()/X_w.cols();

    Eigen::Vector3d T_prenorm(T_prenorm_x,T_prenorm_y,T_prenorm_z);
    Eigen::Matrix4d DP;
    DP.block<3,3>(0,0) = Eigen::DiagonalMatrix<double,3,3>(1,1,1);
    DP.block<3,1>(0,3) = -T_prenorm;
    DP.block<1,4>(3,0) = Eigen::Vector4d(0,0,0,1).transpose();

    Eigen::MatrixXd X_W_T = DP*X_w;

    // Anisotropic scaling: mean distance to origin should be sqrt(3)
    Eigen::MatrixXd X_W_ABS = X_W_T.array().abs().matrix();
    double aniso_x = X_W_ABS.block(0,0,1,X_W_ABS.cols()).array().sum()/X_W_ABS.cols();
    double aniso_y = X_W_ABS.block(1,0,1,X_W_ABS.cols()).array().sum()/X_W_ABS.cols();
    double aniso_z = X_W_ABS.block(2,0,1,X_W_ABS.cols()).array().sum()/X_W_ABS.cols();

    Eigen::Vector3d aniso(aniso_x,aniso_y,aniso_z);
    Eigen::Vector3d ones(1.0,1.0,1.0);

    Eigen::Vector3d prescaleVector = ones.array()/aniso.array();
    Eigen::DiagonalMatrix<double,4> diag;
    diag.diagonal() = Eigen::Vector4d(prescaleVector(0),prescaleVector(1),prescaleVector(2),1);

    Eigen::Matrix4d DS = diag;




    // Combine translation and scaling
    Eigen::Matrix4d DSP = DS*DP;


    // Apply prenormalization

    X_w = DSP*X_w.eval();



    //----------------- Create line measurements vectors-------------------------------------------------------------//

    Eigen::Matrix<double, 3, Eigen::Dynamic> x_1 = Eigen::MatrixXd::Map(x_c.data(),6, x_c.cols()/2).topRows(3);
    Eigen::Matrix<double, 3, Eigen::Dynamic> x_2 = Eigen::MatrixXd::Map(x_c.data(),6, x_c.cols()/2).bottomRows(3);

    Eigen::Matrix<double, 3, Eigen::Dynamic> l_c;
    l_c.resize(3, nlines);

    for(int i = 0; i < nlines; i++){
      l_c.block<3,1>(0,i) = x_1.block<3,1>(0,i).cross(x_2.block<3,1>(0,i));
    }

    //-----------------Prenormalize lines----------------------------------------------------------------------------//



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
    Eigen::MatrixXd M_pl_two = kron(X_1_X_2.transpose(), Eigen::Vector3d(1.0,1.0,1.0).transpose());

    Eigen::MatrixXd M_pl = (M_pl_one.array()*M_pl_two.array()).matrix();


    //------------------------------------ Pose Estimation ------------------------//

    // Perform a singular value decomposition
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M_pl, Eigen::ComputeFullV);
    Eigen::MatrixXd V_right = svd.matrixV().rightCols(1);
    V_right.resize(3,4);

    // First estimate of projection matrix using last singular vector
    Eigen::Matrix<double, 3, 4> P_est = V_right;

    // Revert point prenormalization
    P_est = P_est*DS;

    // Algorithm 1 + 2 from paper: scale estimation + R1 orthogonalization
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_P(P_est.leftCols(3), Eigen::ComputeFullV|Eigen::ComputeFullU);
    double det = (svd_P.matrixU()*((svd_P.matrixV()).transpose())).determinant();
    double scale = 3*det/(svd_P.singularValues().array().sum());


    std::cout << "Condition Factor = " << ((svd_P.singularValues())(0,0))/(svd_P.singularValues()(2,0)) <<"\n\n";


    Eigen::MatrixXd R_est = det*(svd_P.matrixU()*svd_P.matrixV().transpose());
    Eigen::MatrixXd T_est = scale*R_est.transpose()*P_est.rightCols(1)-T_prenorm;

    tf.block<3,3>(0,0) = R_est;
    tf.block<3,1>(0,3) = R_est*T_est;
    tf.block<1,4>(3,0) = Eigen::Vector4d(0,0,0,1).transpose();

    std::cout << "Estimated Pose = \n\n"<<tf<<"\n\n";


  }
  toc();
  return tf;


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
    L.block<3,1>(0,i) = X_1.block<3,1>(0,i).cross(X_2.block<3,1>(0,i));
    L.block<3,1>(3,i) = X_2.block<3,1>(0,i)-X_1.block<3,1>(0,i);
  }

  //std::cout << "3D Line Matrix = \n\n" << L << "\n\n";

  return L;

}


 inline Eigen::Matrix3d Alined::skew(const Eigen::Vector3d& vec){

   Eigen::Matrix3d skew;
   skew << 0 , -vec(2), vec(1),vec(2),0,-vec(0),-vec(1),vec(0),0;

   return skew;
 }

 inline Eigen::Vector3d Alined::unskew(const Eigen::Matrix3d& skew){

   Eigen::Vector3d vec;

   vec << skew(2,1), skew(0,2), skew(1,0);

   return vec;

 }

#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <string>
#include <fstream>


#include <Eigen/Eigen>
#include <Eigen/SVD>


double norm_1(Eigen::Matrix3d X)
{
  double sum = 0;
  for ( int i = 0; i < 3; i++ )
  {
    for ( int j = 0; j < 3; j++ )
    {
      sum += abs(X(i,j));
    }
  }
  return sum;
}

double norm_inf(Eigen::Matrix3d X)
{
  double max_abs = abs(X(0,0));
  for ( int i = 0; i < 3; i++ )
  {
    for ( int j = 0; j < 3; j++ )
    {
      double new_abs = abs(X(i,j));
      if ( new_abs > max_abs )
      {
        max_abs = new_abs;
      }
    }
  }
  return max_abs;
}

template<typename Mat>
void polar_eigen(const Mat& A, Mat& R)
{
  typedef Eigen::Matrix<typename Mat::Scalar,3,1> Vec;
  Eigen::SelfAdjointEigenSolver<Mat> eig;
  eig.computeDirect(A.transpose()*A);
  Vec S = eig.eigenvalues().cwiseSqrt();

  R = A  * eig.eigenvectors() * S.asDiagonal().inverse()
    * eig.eigenvectors().transpose();
}


int main() {
	
  std::ifstream file;
  file.open("SVD_benchmark");
  if (!file) 
  {
    CGAL_TRACE_STREAM << "Error loading file!\n";
    return 0;
  }

  int ite = 200000;
  Eigen::JacobiSVD<Eigen::Matrix3d> svd;
  Eigen::Matrix3d A, U, H;           

  int matrix_idx = rand()%200;
  for (int i = 0; i < matrix_idx; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        file >> A(j, k);
      }
    }
  }

  CGAL::Timer task_timer; 

  CGAL_TRACE_STREAM << "Start SVD decomposition...";
  task_timer.start();
  for (int i = 0; i < ite; i++)
  {
    polar_eigen<Eigen::Matrix3d> (A, U);
    U = U.transpose();
  }
  task_timer.stop();
  file.close();

  CGAL_TRACE_STREAM << "done: " << task_timer.time() << "s\n";

	return 0;
}


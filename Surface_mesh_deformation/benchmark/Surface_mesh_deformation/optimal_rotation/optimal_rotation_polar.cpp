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


void algorithm_polar(Eigen::Matrix3d A, Eigen::Matrix3d &U, double tolerance)
{
  Eigen::Matrix3d X = A;
  int k = -1;
  Eigen::Matrix3d Y;
  double alpha, beta, gamma;
  do 
  {
    k++;
    Y = X.inverse();
    alpha = sqrt( norm_1(X) * norm_inf(X) );
    beta = sqrt( norm_1(Y) * norm_inf(Y) );
    gamma = sqrt(beta/alpha);
    X = 0.5*( gamma*X + Y.transpose()/gamma );

  } while ( abs(gamma-1) > tolerance );

  U = X;
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
  file.close();

  A(0,0) = 0.0014667022958104185; A(0,1) = 0.0024644551180079627; A(0,2) = 0.0040659871060825092;
  A(1,0) = 0.0028327558016991478; A(1,1) = 0.0054236820249146406; A(1,2) = 0.0079280090866983826;
  A(2,0) = 0.0039090073251031232; A(2,1) = 0.0066744523074963374; A(2,2) = 0.010848552426718550;
  A = A*10000;

  CGAL::Timer task_timer; 

  CGAL_TRACE_STREAM << "Start SVD decomposition...";
  task_timer.start();
  for (int i = 0; i < ite; i++)
  {
    algorithm_polar(A, U, 1e-6);
    U.transposeInPlace();
    double det_2 = U.determinant();
    int aaa = 0;
  }
  task_timer.stop();
  

  CGAL_TRACE_STREAM << "done: " << task_timer.time() << "s\n";

	return 0;
}


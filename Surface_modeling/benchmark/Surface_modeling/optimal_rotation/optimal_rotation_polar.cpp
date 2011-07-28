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
    algorithm_polar(A, U, 1e-6);
    U = U.transpose();
  }
  task_timer.stop();
  file.close();

  CGAL_TRACE_STREAM << "done: " << task_timer.time() << "s\n";

	return 0;
}


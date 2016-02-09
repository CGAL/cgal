
#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <string>
#include <fstream>


#include <Eigen/Eigen>
#include <Eigen/SVD>

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
  Eigen::Matrix3d u, v, cov, r;         
  Eigen::Vector3d w;   

  int matrix_idx = rand()%200;
  for (int i = 0; i < matrix_idx; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        file >> cov(j, k);
      }
    }
  }


  CGAL::Timer task_timer; 

  CGAL_TRACE_STREAM << "Start SVD decomposition...";
  task_timer.start();
  for (int i = 0; i < ite; i++)
  {
    
    svd.compute( cov, Eigen::ComputeFullU | Eigen::ComputeFullV );
    u = svd.matrixU(); v = svd.matrixV(); w = svd.singularValues();
    r = v*u.transpose();
  }
  task_timer.stop();
  file.close();

  CGAL_TRACE_STREAM << "done: " << task_timer.time() << "s\n";

  return 0;
}
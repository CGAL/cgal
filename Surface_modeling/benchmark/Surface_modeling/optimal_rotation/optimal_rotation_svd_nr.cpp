#include "svd.h"
#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <string>
#include <fstream>

int main() {
	
  std::ifstream file;
  file.open("SVD_benchmark");
  if (!file) 
  {
    CGAL_TRACE_STREAM << "Error loading file!\n";
    return 0;
  }

  // initialization of matrices
  int ite = 200000;
  double ***u;
  u = new double** [ite];
  for (int i = 0; i < ite; i++)
  {
    u[i] = new double* [4];
    for (int j = 0; j < 4; j++)
    {
      u[i][j] = new double [4];
    }
  }
  int matrix_idx = rand()%200;
  for (int i = 0; i < matrix_idx; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        file >> u[0][j+1][k+1];
      }
    }
  }
  for (int i = 1; i < ite; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        u[i][j+1][k+1] = u[0][j+1][k+1];
      }
    }
  }
  double **v;
  v = new double* [4];
  for (int i = 0; i < 4; i++)
  {
    v[i] = new double[4];
  }
  double *w;
  w = new double[4];

  CGAL::Timer task_timer; 

  CGAL_TRACE_STREAM << "Start SVD decomposition...";
  task_timer.start();
  for (int i = 0; i < ite; i++)
  {
    svdcmp(u[i], 3, 3, w, v);
  }
  task_timer.stop();
  file.close();

  CGAL_TRACE_STREAM << "done: " << task_timer.time() << "s\n";

  delete [] w;
  for (int i = 0; i < 4; i++)
  {
    delete [] v[i];
  }
  delete [] v;
  for (int i = 0; i < ite; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      delete [] u[i][j];
    }
    delete [] u[i];
  }
  delete [] u;

	return 0;
}


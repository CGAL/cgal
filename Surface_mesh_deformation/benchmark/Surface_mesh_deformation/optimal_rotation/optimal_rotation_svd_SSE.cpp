#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <string>
#include <fstream>


#include <Eigen/Eigen>

// #define USE_SCALAR_IMPLEMENTATION
#define USE_SSE_IMPLEMENTATION
#define COMPUTE_V_AS_MATRIX
#define COMPUTE_U_AS_MATRIX

#include <CGAL/internal/Surface_mesh_deformation/auxiliary/Singular_Value_Decomposition_Preamble.hpp>

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
  Eigen::Matrix3d u, v, m, r;
  Eigen::Vector3d w;

  int matrix_idx = rand()%200;
  for (int i = 0; i < matrix_idx; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        file >> m(j, k);
      }
    }
  }


  CGAL::Timer task_timer;

  std::cout << "Start SVD decomposition...";
  task_timer.start();
  for (int i = 0; i < ite; i++)
  {

    #include <CGAL/internal/Surface_mesh_deformation/auxiliary/Singular_Value_Decomposition_Kernel_Declarations.hpp>

    ENABLE_SCALAR_IMPLEMENTATION(Sa11.f=m(0,0);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa21.f=m(1,0);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa31.f=m(2,0);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa12.f=m(0,1);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa22.f=m(1,1);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa32.f=m(2,1);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa13.f=m(0,2);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa23.f=m(1,2);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa33.f=m(2,2);)

    ENABLE_SSE_IMPLEMENTATION(Va11=_mm_set1_ps(m(0,0));)
    ENABLE_SSE_IMPLEMENTATION(Va21=_mm_set1_ps(m(1,0));)
    ENABLE_SSE_IMPLEMENTATION(Va31=_mm_set1_ps(m(2,0));)
    ENABLE_SSE_IMPLEMENTATION(Va12=_mm_set1_ps(m(0,1));)
    ENABLE_SSE_IMPLEMENTATION(Va22=_mm_set1_ps(m(1,1));)
    ENABLE_SSE_IMPLEMENTATION(Va32=_mm_set1_ps(m(2,1));)
    ENABLE_SSE_IMPLEMENTATION(Va13=_mm_set1_ps(m(0,2));)
    ENABLE_SSE_IMPLEMENTATION(Va23=_mm_set1_ps(m(1,2));)
    ENABLE_SSE_IMPLEMENTATION(Va33=_mm_set1_ps(m(2,2));)

    #include <CGAL/internal/Surface_mesh_deformation/auxiliary/Singular_Value_Decomposition_Main_Kernel_Body.hpp>

    std::pair<Eigen::Matrix3d, Eigen::Matrix3d> solver;

#ifdef USE_SCALAR_IMPLEMENTATION
    solver.first(0,0) = Su11.f;
    solver.first(1,0) = Su21.f;
    solver.first(2,0) = Su31.f;
    solver.first(0,1) = Su12.f;
    solver.first(1,1) = Su22.f;
    solver.first(2,1) = Su32.f;
    solver.first(0,2) = Su13.f;
    solver.first(1,2) = Su23.f;
    solver.first(2,2) = Su33.f;

    solver.second(0,0) = Sv11.f;
    solver.second(1,0) = Sv21.f;
    solver.second(2,0) = Sv31.f;
    solver.second(0,1) = Sv12.f;
    solver.second(1,1) = Sv22.f;
    solver.second(2,1) = Sv32.f;
    solver.second(0,2) = Sv13.f;
    solver.second(1,2) = Sv23.f;
    solver.second(2,2) = Sv33.f;
#endif

#ifdef USE_SSE_IMPLEMENTATION
    float buf[4];
    _mm_storeu_ps(buf,Vu11);solver.first(0,0)=buf[0];
    _mm_storeu_ps(buf,Vu21);solver.first(1,0)=buf[0];
    _mm_storeu_ps(buf,Vu31);solver.first(2,0)=buf[0];
    _mm_storeu_ps(buf,Vu12);solver.first(0,1)=buf[0];
    _mm_storeu_ps(buf,Vu22);solver.first(1,1)=buf[0];
    _mm_storeu_ps(buf,Vu32);solver.first(2,1)=buf[0];
    _mm_storeu_ps(buf,Vu13);solver.first(0,2)=buf[0];
    _mm_storeu_ps(buf,Vu23);solver.first(1,2)=buf[0];
    _mm_storeu_ps(buf,Vu33);solver.first(2,2)=buf[0];

    _mm_storeu_ps(buf,Vv11);solver.second(0,0)=buf[0];
    _mm_storeu_ps(buf,Vv21);solver.second(1,0)=buf[0];
    _mm_storeu_ps(buf,Vv31);solver.second(2,0)=buf[0];
    _mm_storeu_ps(buf,Vv12);solver.second(0,1)=buf[0];
    _mm_storeu_ps(buf,Vv22);solver.second(1,1)=buf[0];
    _mm_storeu_ps(buf,Vv32);solver.second(2,1)=buf[0];
    _mm_storeu_ps(buf,Vv13);solver.second(0,2)=buf[0];
    _mm_storeu_ps(buf,Vv23);solver.second(1,2)=buf[0];
    _mm_storeu_ps(buf,Vv33);solver.second(2,2)=buf[0];
#endif

    r = solver.first * solver.second.transpose();
  }
  task_timer.stop();
  file.close();

  std::cout << "done: " << task_timer.time() << "s\n";
  return 0;
}
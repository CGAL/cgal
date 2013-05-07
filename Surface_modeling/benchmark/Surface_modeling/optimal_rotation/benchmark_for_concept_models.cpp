#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>
#include <vector>

// models of DeformationClosestRotationTraits_3
#include <CGAL/Deformation_Eigen_closest_rotation_traits_3.h>
#include <CGAL/internal/Surface_modeling/Deformation_Eigen_polar_closest_rotation_traits_3.h>
#include <CGAL/internal/Surface_modeling/auxiliary/Deformation_fast_SVD_closest_rotation_traits_3.h>

// provide a model of DeformationClosestRotationTraits_3
//
// Warning: we require more functionality than requested in DeformationClosestRotationTraits_3,
// so this benchmark may not compile with other models of DeformationClosestRotationTraits_3
template<class DeformationClosestRotationTraits_3>
void benchmark()
{
  typedef typename DeformationClosestRotationTraits_3::Matrix Matrix;

  DeformationClosestRotationTraits_3 model;
  Matrix m = model.zero_matrix(); // no default const is required in concept (this is a workaround)
  Matrix R = model.zero_matrix();

  std::ifstream file("SVD_benchmark");
  if (!file) {
    std::cerr << "Error: can not open SVD_benchmark file!" << std::endl;
    return;
  }

  std::vector<Matrix> matrices;

  while( true ) {    
    bool read_ok = true;
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        if( !(file >> m(j, k)) ) // Warning: this part violates concept (in concept we do not have access to matrix coeffs)// 
        { read_ok = false; break; }       
      }
    }

    if(!read_ok) { break; }
    matrices.push_back(m);
  }
  std::cerr << "Reading matrices from file is completed ( "<< matrices.size() << " number of matrices )." << std::endl;

  std::cerr << "Starting benchmark for " << matrices.size() << " matrices..." << std::endl;
  CGAL::Timer task_timer; task_timer.start();
  for(typename std::vector<Matrix>::iterator it = matrices.begin();
    it != matrices.end(); ++it) {
      model.compute_close_rotation(*it, R);
  }
  std::cerr << "-----------------------------" << std::endl;
  std::cerr << "Done: " << task_timer.time() << " sc" << std::endl;
  std::cerr << "-----------------------------" << std::endl;
}


int main() {

  std::cerr << "Benchmark for Deformation_Eigen_closest_rotation_traits_3 (Eigen SVD): " << std::endl;
  benchmark<CGAL::Deformation_Eigen_closest_rotation_traits_3>();

  std::cerr << "Benchmark for Deformation_Eigen_polar_closest_rotation_traits_3 (Eigen polar(as a filter) and SVD(as a gold standard)): " << std::endl;
  benchmark<CGAL::Deformation_Eigen_polar_closest_rotation_traits_3>();

  std::cerr << "Benchmark for Deformation_fast_SVD_closest_rotation_traits_3 (fast SVD code with SSE): " << std::endl;
  benchmark<CGAL::Deformation_fast_SVD_closest_rotation_traits_3>();

  std::cerr << "All done!" << std::endl;
}
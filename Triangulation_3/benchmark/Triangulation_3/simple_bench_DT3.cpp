//#define CGAL_PROFILE
//#define CGAL_USE_SSE2_FABS
//#define CGAL_USE_SSE2_MAX
//#define CGAL_MSVC_USE_STD_FABS  // use this one with precise
//#define INDEX_STORAGE 1

#include "CGAL/Bbox_3.h"
#include "CGAL/Triangulation_data_structure_3.h"
#include "CGAL/Triangulation_vertex_base_3.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/Memory_sizer.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <string>
#include <fstream>
#include <clocale>

#if PARALLEL
using Concurrent_tag = CGAL::Parallel_tag;
#else
using Concurrent_tag = CGAL::Sequential_tag;
#endif
typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
#if INDEX_STORAGE
typedef CGAL::VertexWithPoint<K> Vb;
typedef CGAL::Cell4Delaunay<K> Cb;
typedef CGAL::Index_tag Tds_type_tag;
#else
typedef CGAL::Triangulation_vertex_base_3<K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb;
typedef CGAL::Handle_tag Tds_type_tag;
#endif
typedef CGAL::Triangulation_data_structure_3<Vb, Cb, Concurrent_tag, Tds_type_tag> Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds>                DT;
typedef DT::Point                                            Point_3;
typedef CGAL::Timer                                          Timer;
typedef CGAL::Memory_sizer                                   Memory_sizer;

int main(int argc, char* argv[])
{
  std::locale loc = std::locale()
      .combine<std::numpunct<char>>(std::locale("en_US.UTF8"));
  std::cout.imbue(loc);
  Memory_sizer memory_sizer;
  auto res_mem_at_start = memory_sizer.resident_size();
  std::cout << "Memory usage right at start:\n" << memory_sizer.virtual_size() << " bytes (virtual), "
            << res_mem_at_start << " bytes (resident)" << std::endl;

  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/ocean_r.xyz");
  std::ifstream in(filename.c_str());
  std::vector<Point_3> points;
  Point_3 p, q;

  while(in >> p ){
    points.push_back(p);
  }

  std::cout << points.size() << " points read\n";

  Timer timer;
  timer.start();
  int M = 100; // Number of times to compute the triangulation
  std::cout << "Compute triangulation "<< M << " times" << std::endl;
  for(int i = 0; i < M; i++){
#if PARALLEL
    DT::Lock_data_structure lock_data_structure{CGAL::bbox_3(points.begin(), points.end()), 50};
    DT dt(&lock_data_structure);
#else
    DT dt;
#endif
    dt.insert(points.begin(), points.end());
    std::cout << "Number of cells: " << dt.number_of_cells() << std::endl;
    auto res_mem = memory_sizer.resident_size();
    std::cout << "Memory usage after construction of the triangulation:\n" << memory_sizer.virtual_size() << " bytes (virtual), "
            << res_mem << " bytes (resident)" << std::endl;
    std::cout << "Diff in resident memory: "
              << res_mem - res_mem_at_start << " bytes" << std::endl;
  }

  timer.stop();
  std::cout << "Time elapsed: " << timer.time() << " sec" << std::endl;

  auto res_mem = memory_sizer.resident_size();
  std::cout << "Memory usage after deallocation of the triangulation:\n" << memory_sizer.virtual_size() << " bytes (virtual), "
            << res_mem << " bytes (resident)" << std::endl;
  std::cout << "Diff in resident memory: "
            << res_mem - res_mem_at_start << " bytes" << std::endl;

  return 0;
}

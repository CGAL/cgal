#include <CGAL/approximate_convex_segmentation.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Real_timer.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> EPICK_Mesh;
typedef boost::graph_traits<EPICK_Mesh>::face_descriptor face_descriptor;

typedef CGAL::Real_timer Timer;

#ifndef CGAL_LINKED_WITH_TBB
typedef CGAL::Sequential_tag Concurrency_tag;
#else
typedef CGAL::Parallel_tag Concurrency_tag;
#endif

template<class Mesh>
double compute_and_time_concavity(const Mesh& mesh)
{
  Timer timer;
  
  timer.start();
  double concavity = CGAL::concavity_values<Concurrency_tag>(mesh);
  timer.stop();

  std::cout << "  concavity value: " << concavity << std::endl;
  std::cout << "  Concavity Time: " << timer.time() << " seconds" << std::endl;
  return timer.time();
}

template <class Mesh>
double compute_and_time_segmentation(const Mesh& mesh, double concavity_threshold)
{
  typedef std::map<face_descriptor, std::size_t> Segments_id_map;
  Segments_id_map segments_map;
  typedef boost::associative_property_map<Segments_id_map> Segments_id_pmap;
  Segments_id_pmap segments_pmap(segments_map);

  Timer timer;
  timer.start();
  std::size_t number_of_segments = CGAL::approximate_convex_segmentation<Concurrency_tag>(mesh,
                                                                                          segments_pmap,
                                                                                          concavity_threshold);
  timer.stop();
  std::cout << "  number of segments: " << number_of_segments << std::endl;
  std::cout << "  Segmentation Time with concavity_threshold=" << concavity_threshold << " : " << timer.time() << " seconds" << std::endl;
  return timer.time();
}

template <class Mesh>
std::vector<double> read_and_run(const std::string& file_name)
{
  std::cout << "File name " << file_name;
  std::vector<double> timings;
  
  Mesh mesh;
  std::ifstream input(file_name.c_str());
  if (!input || !(input >> mesh) || CGAL::is_empty(mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Not a valid off file." << std::endl;
    assert(false);
  }
  std::cout << " Number of triangles: " << CGAL::num_faces(mesh) << std::endl;

  timings.push_back(CGAL::num_faces(mesh)); // [0]

  double concavity_time = compute_and_time_concavity(mesh);
  timings.push_back(concavity_time);

  std::vector<double> concavity_threshold_values;
  concavity_threshold_values.push_back(0.01);
  concavity_threshold_values.push_back(0.05);
  concavity_threshold_values.push_back(0.5);

  for (std::size_t i = 0; i < concavity_threshold_values.size(); ++i)
  {
    double segment_time = compute_and_time_segmentation(mesh, concavity_threshold_values[i]);
    timings.push_back(segment_time); // [1..]
  }

  return timings;
}

int main()
{
  std::vector<std::string> files;
  files.push_back("data/sword.off");
  files.push_back("data/cactus.off");
  files.push_back("data/elephant.off");
  files.push_back("data/dino.off");

  std::vector<std::vector<double> > timings_epick;

  for (std::size_t i = 0; i < files.size(); ++i)
  {
    timings_epick.push_back(read_and_run<EPICK_Mesh>(files[i]));
  }

  // print timings
  std::cout << std::endl << std::endl;
  std::cout << "*** Concavity values | Timings (EPICK) ****" << std::endl;
  std::cout << "+------------+---------------+" << std::endl;
  std::cout << "| #Triangles | Run-time (ms) |" << std::endl;
  std::cout << "+------------+---------------+" << std::endl;
  for (std::size_t i = 0; i < timings_epick.size(); ++i)
  {
    std::cout << "| "
              << std::setw(11) << timings_epick[i][0] << "|"
              << std::setw(15)  << timings_epick[i][1] * 1000.  << "|" << std::endl; 
  }
  std::cout << "+------------+---------------+" << std::endl;

  std::cout << std::endl << std::endl;
  std::cout << "*** Approximate convex segmentation | Timings (EPICK) ****" << std::endl;
  std::cout << "+------------+-----------------+-----------------+----------------+" << std::endl;
  std::cout << "| #Triangles | #Concavity=0.01 | #Concavity=0.05 | #Concavity=0.5 |" << std::endl;
  std::cout << "+------------+-----------------+-----------------+----------------+" << std::endl;
  for (std::size_t i = 0; i < timings_epick.size(); ++i)
  {
    std::cout << "| " 
              << std::setw(11) << timings_epick[i][0] << "|"
              << std::setw(17) << timings_epick[i][2] << "|"
              << std::setw(17) << timings_epick[i][3] << "|"
              << std::setw(16) << timings_epick[i][4] << "|" << std::endl;
  }
  std::cout << "+------------+-----------------+-----------------+----------------+" << std::endl;
}

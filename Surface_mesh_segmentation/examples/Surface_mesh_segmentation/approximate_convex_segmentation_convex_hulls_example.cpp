#include <CGAL/approximate_convex_segmentation.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Real_timer.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

typedef CGAL::Real_timer Timer;

#ifndef CGAL_LINKED_WITH_TBB
typedef CGAL::Sequential_tag Concurrency_tag;
#else
typedef CGAL::Parallel_tag Concurrency_tag;
#endif

int main()
{
  // read mesh
  Mesh mesh;
  
  std::ifstream input("data/sword.off");
  
  if (!input || !(input >> mesh))
  {
    std::cout << "Failed to read mesh" << std::endl;
    return EXIT_FAILURE;
  }

  if (CGAL::is_empty(mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cout << "Input mesh is invalid" << std::endl;
    return EXIT_FAILURE;
  }

  // create property map for segment-ids
  typedef Mesh::Property_map<face_descriptor, int> Clusters_pmap;
  Clusters_pmap segments_pmap = mesh.add_property_map<face_descriptor, int>("f:segment").first;

  // create property map for convex hulls
  boost::vector_property_map<Mesh> convex_hulls_pmap;

  // decompose mesh
  Timer timer;
  
  timer.start();
  std::size_t segments_num = CGAL::approximate_convex_segmentation<Concurrency_tag>(mesh, segments_pmap, 0.3, CGAL::parameters::segments_convex_hulls(convex_hulls_pmap));
  timer.stop();

  std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;

  // write segment-ids for each facet
  std::cout << "Number of segments: " << segments_num << std::endl;
  BOOST_FOREACH(face_descriptor fd, faces(mesh))
  {
    std::cout << segments_pmap[fd] << " ";
  }
  std::cout << std::endl;

  // write concavity values for all segments
  for (std::size_t i = 0; i < segments_num; ++i)
  {
    std::cout << "Concavity value of #" << i << " segment: " << CGAL::concavity_values<Concurrency_tag>(mesh, segments_pmap, i) << std::endl;
  }

  // use convex hulls
  for (std::size_t i = 0; i < segments_num; ++i)
  {
    Mesh& convex_hull = convex_hulls_pmap[i];
    // ...
  }

  return EXIT_SUCCESS;
}


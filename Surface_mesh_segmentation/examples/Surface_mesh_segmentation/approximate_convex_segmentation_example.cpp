#include <CGAL/approximate_convex_segmentation.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Real_timer.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;
typedef typename boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

typedef CGAL::Real_timer Timer;

#ifndef CGAL_LINKED_WITH_TBB
typedef CGAL::Sequential_tag Concurrency_tag;
#else
typedef CGAL::Parallel_tag Concurrency_tag;
#endif

int main()
{
  // read mesh
  Polyhedron mesh;

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
 
  // init the polyhedron indices
  CGAL::set_halfedgeds_items_id(mesh);

  // create property map for segment-ids
  typedef boost::property_map<Polyhedron, CGAL::dynamic_face_property_t<int> >::type Facet_property_map;
  Facet_property_map facet_property_map = get(CGAL::dynamic_face_property_t<int>(), mesh);

  // decompose mesh with default parameters
  Timer timer;

  timer.start(); 
  std::size_t segments_num = CGAL::approximate_convex_segmentation<Concurrency_tag>(mesh, facet_property_map, 0.3);
  timer.stop();

  std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;

  // write segment-ids for each facet
  std::cout << "Number of segments: " << segments_num << std::endl;
  BOOST_FOREACH(face_descriptor face, faces(mesh))
  {
    std::cout << get(facet_property_map, face) << " ";
  }
  std::cout << std::endl;

  // write concavity values for all segments
  for (std::size_t i = 0; i < segments_num; ++i)
  {
    std::cout << "Concavity value of #" << i << " segment: " << CGAL::concavity_values<Concurrency_tag>(mesh, facet_property_map, i) << std::endl;
  }

  return EXIT_SUCCESS;
}


#include <CGAL/approximate_convex_segmentation.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "Utils.h" 

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

#ifndef CGAL_LINKED_WITH_TBB
typedef CGAL::Sequential_tag Concurrency_tag;
#else
typedef CGAL::Parallel_tag Concurrency_tag;
#endif

int main()
{
  Mesh mesh;
  if (!read_to_polyhedron("data/sword.off", mesh)) return EXIT_FAILURE;

  typedef Mesh::Property_map<face_descriptor, int> Clusters_map;
  Clusters_map segment_ids = mesh.add_property_map<face_descriptor, int>("f:segment").first;

  std::size_t segments_num = CGAL::approximate_convex_segmentation<Concurrency_tag>(mesh, segment_ids, 0.3);

  std::cout << "Number of segments: " << segments_num << std::endl;
  BOOST_FOREACH(face_descriptor fd, faces(mesh))
  {
    std::cout << segment_ids[fd] << " ";
  }
  std::cout << std::endl;

  for (std::size_t i = 0; i < segments_num; ++i)
  {
    double concavity = CGAL::concavity_values<Concurrency_tag>(mesh, segment_ids, i);
    std::cout << "Concavity value of #" << i << " segment: " << concavity << std::endl;

    if (concavity < 0 || concavity > 0.3)
    {
      std::cerr << "Resulting concavity value doesn't meet constraints" << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}


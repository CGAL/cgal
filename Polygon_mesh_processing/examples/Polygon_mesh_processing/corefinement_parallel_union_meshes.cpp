#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  std::ifstream in(argc>1?argv[1]:"data/eight.off");
  int nb_copies = argc > 2 ? atoi(argv[2]) : 100;
  if (nb_copies<=0) return 1;
  std::vector<Mesh> meshes(nb_copies);
  in >> meshes[0];
  for (int i=1; i<nb_copies; ++i)
  {
    CGAL::Aff_transformation_3<Kernel> trans(CGAL::Translation(),  Kernel::Vector_3(i*0.2, 0, 0));
    meshes[i]=meshes[0];
    PMP::transform(trans, meshes[i]);
  }

  // lambda function used to do the union of two meshes:
  // we join the left part of the mesh vector with the right part
  // as the right part will be discarded by the resize in the while loop
  auto f = [&meshes](const tbb::blocked_range<std::size_t>& range)
           {
             for( std::size_t k = range.begin(); k != range.end(); ++k)
             {
               PMP::corefine_and_compute_union(meshes[k], meshes[meshes.size() - k - 1], meshes[k]);
             }
           };

  // do the union of meshes in parallel
  while (meshes.size()>1)
  {
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, meshes.size()/2), f);
    meshes.resize((meshes.size() + 1) / 2);
  }

  std::ofstream out("multiple_union.off");
  out << std::setprecision(17) << meshes[0];

  return 0;
}

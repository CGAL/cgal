#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Polyhedron_3<Kernel> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type Vertex_point_pmap;

typedef CGAL::Search_traits_3<Kernel>                                    Traits_base;
typedef CGAL::Search_traits_adapter<vertex_descriptor, Vertex_point_pmap, Traits_base> Traits;

typedef CGAL::Kd_tree<Traits>                      Tree;
typedef Tree::Splitter                             Splitter;

int main(int argc, char* argv[])
{
  Mesh mesh;
  std::ifstream in((argc>1)?argv[1]:"data/tripod.off");
  in  >> mesh;

  Vertex_point_pmap vppmap = get(CGAL::vertex_point,mesh);
  Traits traits(vppmap);
  Tree tree(vertices(mesh).begin(), vertices(mesh).end(), Splitter(), traits);

  Point_3 query(0.0, 0.0, 0.0);
  double radius = 0.5;
  double epsilon = 0.01;

  // search vertices
  CGAL::Fuzzy_sphere<Traits> fz(query, radius, epsilon, traits);

  //collect vertices that are inside the sphere
  std::list<vertex_descriptor> result;
  tree.search(std::back_inserter(result), fz);
  std::cout << "There are " << result.size() << " vertices inside the fuzzy sphere\n";

  return 0;
}

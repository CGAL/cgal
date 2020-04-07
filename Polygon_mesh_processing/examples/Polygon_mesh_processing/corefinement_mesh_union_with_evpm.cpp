#define CGAL_PROFILE
#define CGAL_USE_FILTERED_RATIONAL_KERNEL
#define MSC_USE_DLL 1
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Timer.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_kernel;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

struct Exact_vertex_point_map
{
  typedef Mesh::Property_map<Mesh::Vertex_index, Exact_kernel::Point_3> Exact_vpm;
  typedef Mesh::Property_map<Mesh::Vertex_index, bool> Exact_vpm_initialized;
  typedef Mesh::Vertex_index key_type;
  typedef Exact_kernel::Point_3 value_type;
  typedef const value_type& reference;
  typedef boost::lvalue_property_map_tag category;

  Exact_vertex_point_map()
    : tm_ptr(nullptr)
  {}


  Exact_vertex_point_map(Mesh& tm)
    : tm_ptr(&tm)
  {
    bool is_new=false;
    std::tie(e_vpm, is_new) =
      tm.add_property_map<Mesh::Vertex_index, Exact_kernel::Point_3>("v:exact_vpm");
    if (is_new)
    {
      CGAL::Cartesian_converter<K, Exact_kernel> to_exact;
      for (Mesh::Vertex_index v : tm_ptr->vertices())
        e_vpm[v] = to_exact( tm_ptr->point(v) );
    }
  }

  friend
  reference
  get(const Exact_vertex_point_map& m, Mesh::Vertex_index v)
  {
    CGAL_assertion(m.tm_ptr!=NULL);
    return m.e_vpm[v];
  }

  friend
  void
  put(const Exact_vertex_point_map& m, Mesh::Vertex_index v, const Exact_kernel::Point_3& p)
  {
    CGAL_assertion(m.tm_ptr!=NULL);
    m.e_vpm[v] = p;
    m.tm_ptr->point(v) = m.to_input( p );
  }

  Mesh* tm_ptr;
  Exact_vpm e_vpm;
  CGAL::Cartesian_converter<Exact_kernel, K> to_input;
};

int main(int argc, char* argv[])
{
  CGAL::Timer t;
  {
  const char* filename1 = (argc > 1) ? argv[1] : "data/blobby.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/eight.off";
  std::ifstream input(filename1);

  Mesh mesh1, mesh2;
  if (!input || !(input >> mesh1))
  {
    std::cerr << "First mesh is not a valid off file." << std::endl;
    return 1;
  }
  input.close();
  input.open(filename2);
  if (!input || !(input >> mesh2))
  {
    std::cerr << "Second mesh is not a valid off file." << std::endl;
    return 1;
  }


  t.start();
  Mesh out;
  bool valid_union = PMP::corefine_and_compute_union(mesh1, mesh2, out,
                                                     CGAL::parameters::vertex_point_map(Exact_vertex_point_map(mesh1)),
                                                     CGAL::parameters::vertex_point_map(Exact_vertex_point_map(mesh2)),
                                                     CGAL::parameters::vertex_point_map(Exact_vertex_point_map(out)) );


  if (! valid_union)
  {
    std::cout << "Union could not be computed\n";
    return 1;
  }
  std::cout << "Union was successfully computed\n";
  //std::ofstream output("union.off");
  //output << out;

  }
  std::cout << t.time() << " sec." << std::endl;
  return 0;
}

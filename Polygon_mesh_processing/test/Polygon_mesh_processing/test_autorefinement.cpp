#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

template <class TriangleMesh>
struct My_visitor :
  public CGAL::Polygon_mesh_processing::Corefinement::Default_visitor<TriangleMesh>
{
  void after_subface_creations(TriangleMesh&){++(*i);}

  My_visitor()
    : i (new int(0) )
  {}

  std::shared_ptr<int> i;
};

void test(const char* fname, std::size_t nb_polylines, std::size_t total_nb_points,
          std::size_t nb_vertices_after_autorefine, bool all_fixed, std::size_t nb_vertices_after_fix,
          bool triple_intersection)
{
  std::cout << "Running tests on " << fname << "\n";
  std::ifstream input(fname);

  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "  Input mesh is not a valid off file." << std::endl;
    exit(EXIT_FAILURE);
  }
  input.close();
  std::size_t nb_vertices_before_autorefine = num_vertices(mesh);

// Testing surface_self_intersection()
  try{
    std::vector< std::vector<K::Point_3> >polylines;
    PMP::experimental::surface_self_intersection(mesh, std::back_inserter(polylines));
    assert(polylines.size() == nb_polylines);
    std::size_t total_nb_pt=0;
    for(const std::vector<K::Point_3>& polyline : polylines)
      total_nb_pt+=polyline.size();
    assert(total_nb_points == total_nb_pt);
    assert( !triple_intersection );
  }
  catch(const PMP::Corefinement::Triple_intersection_exception&)
  {
    assert( triple_intersection );
  }

// Testing autorefine()
  try{
    My_visitor<Mesh> visitor;
    PMP::experimental::autorefine(mesh,
      PMP::parameters::visitor(visitor));
    mesh.collect_garbage();
    assert( nb_vertices_after_autorefine==num_vertices(mesh));
    assert( (nb_vertices_before_autorefine!=nb_vertices_after_autorefine)== (*(visitor.i) != 0) );
    assert( !triple_intersection );
  }
  catch(const PMP::Corefinement::Triple_intersection_exception&)
  {
    assert( triple_intersection );
  }

// Testing autorefine_and_remove_self_intersections()
  try{
    input.open(fname);
    mesh.clear();
    input >> mesh;
    bool res=PMP::experimental::autorefine_and_remove_self_intersections(mesh);
    assert(res==all_fixed);
    mesh.collect_garbage();
    assert( nb_vertices_after_fix==num_vertices(mesh));
    assert( !triple_intersection );
  }
  catch(const PMP::Corefinement::Triple_intersection_exception&)
  {
    assert( triple_intersection );
  }
}

int main(int argc, const char** argv)
{
  // file nb_polylines total_nb_points nb_vertices_after_autorefine all_fixed nb_vertices_after_fix triple_intersection
  for (int i=0;i<(argc-1)/7; ++i)
    test(argv[1+7*i], atoi(argv[1+7*i+1]), atoi(argv[1+7*i+2]),
         atoi(argv[1+7*i+3]), atoi(argv[1+7*i+4])==0?false:true, atoi(argv[1+7*i+5]), atoi(argv[1+7*i+6])==0?false:true);
}

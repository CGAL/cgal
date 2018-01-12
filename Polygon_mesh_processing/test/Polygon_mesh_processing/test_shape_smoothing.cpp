#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/shape_smoothing.h>
#include <boost/graph/graph_traits.hpp>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

struct Constraints_pmap
{
  std::set<vertex_descriptor>* set_ptr_;

  typedef vertex_descriptor                   key_type;
  typedef bool                                value_type;
  typedef value_type&                         reference;
  typedef boost::read_write_property_map_tag  category;

public:
  Constraints_pmap(std::set<vertex_descriptor>* set_ptr)
    : set_ptr_(set_ptr)
  {}
  Constraints_pmap()
    : set_ptr_(NULL)
  {}

  friend value_type get(const Constraints_pmap& map, const key_type& e)
  {
    CGAL_assertion(map.set_ptr_ != NULL);
    return !map.set_ptr_->empty()
         && map.set_ptr_->count(e);
  }
  friend void put(Constraints_pmap& map
                , const key_type& e, const value_type is)
  {
    CGAL_assertion(map.set_ptr_ != NULL);
    if (is)                map.set_ptr_->insert(e);
    else if(get(map, e))   map.set_ptr_->erase(e);
  }
};

void test_implicit_constrained(const char* filename)
{
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;
  input.close();

  boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  // z max is 20 in the devil;
  std::set<vertex_descriptor> selected_vertices;

  for (vertex_descriptor v : vertices(mesh))
  {
    const double z = get(vpmap, v).z();
    if(z  > 19)
      selected_vertices.insert(v);
  }

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::cout<<"number of selected vertices= "<<selected_vertices.size()<<std::endl;
  #endif

  Constraints_pmap vcmap(&selected_vertices);

  const double time_step = 1;

  CGAL::Polygon_mesh_processing::smooth_modified_curvature_flow(mesh, time_step,
                            CGAL::Polygon_mesh_processing::parameters::vertex_is_constrained_map(vcmap).
                                                                       number_of_iterations(5));
#define CGAL_PMP_SMOOTHING_VERBOSE
  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("data/output_implicit.off");
  out << mesh;
  out.close();
  #endif
}

void test_explicit_scheme(const char* filename)
{
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;
  input.close();

  boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  const double time_step = 1;
  const unsigned int iterations = 5;
  CGAL::Polygon_mesh_processing::smooth_modified_curvature_flow(mesh, time_step,
                            CGAL::Polygon_mesh_processing::parameters::use_explicit_scheme(true).
                                                                       number_of_iterations(iterations));


  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  std::ofstream out("data/output_explicit.off");
  out << mesh;
  out.close();
  #endif
}

int main(int argc, char* argv[])
{
  const char* filename = "data/mannequin-devil.off";
  std::ifstream input(filename);

  test_implicit_constrained(filename);
  test_explicit_scheme(filename);

  return 0;
}

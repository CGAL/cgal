#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Real_timer.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

// BGL bug-fix
namespace CGAL{
namespace Polygon_mesh_processing{
template <typename TriangleMesh, typename NamedParameters>
bool is_degenerate_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                 const TriangleMesh& tm,
                                 const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle(halfedge(f, tm), tm));

  using boost::get_param;
  using boost::choose_param;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, tm));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Traits;
  Traits traits = choose_param(get_param(np, internal_np::geom_traits), Traits());

  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h = halfedge(f, tm);

  return traits.collinear_3_object()(get(vpmap, source(h, tm)),
                                     get(vpmap, target(h, tm)),
                                     get(vpmap, target(next(h, tm), tm)));
}

template <typename TriangleMesh>
bool is_degenerate_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                 const TriangleMesh& tm)
{
  return CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(f, tm, parameters::all_default());
}
}

}

#include <CGAL/boost/graph/io.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// Polyhedron type
typedef CGAL::Surface_mesh<K::Point_3> Polyhedron;
// Domain 
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;


// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,CGAL::Parallel_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int, char** argv)
{
  // Load OBJ
  std::vector<K::Point_3> points;
  std::vector<std::vector<std::size_t> > faces;

  std::ifstream in(argv[1]);
  if(!in || !CGAL::read_OBJ(in,points,faces))
  {
      return 1;
  }

  Polyhedron poly;
  namespace PMP = CGAL::Polygon_mesh_processing;
  PMP::orient_polygon_soup(points,faces);
  PMP::polygon_soup_to_polygon_mesh(points, faces, poly);
  if (!CGAL::is_triangle_mesh(poly))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  if (!CGAL::is_triangle_mesh(poly)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // remove degenerate faces
  std::vector<Polyhedron::Face_index> faces_to_remove;
  for (Polyhedron::Face_index f : poly.faces())
    if( PMP::is_degenerate_triangle_face(f, poly))
      faces_to_remove.push_back(f);
  for (Polyhedron::Face_index f : faces_to_remove)
    CGAL::Euler::remove_face(halfedge(f,poly), poly);
  faces_to_remove.clear();

  // Create a vector with only one element: the pointer to the polyhedron.
  std::vector<Polyhedron*> poly_ptrs_vector(1, &poly);

  // Create a polyhedral domain, with only one polyhedron,
  // and no "bounding polyhedron", so the volumetric part of the domain will be
  // empty.
  Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());
  
  // Get sharp features
  //  domain.detect_features(); //includes detection of borders

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 1,
                         facet_angle = 25,
                         facet_size = 0.2,
                         facet_distance = 0.02);
  CGAL::Real_timer timer;
  timer.start();
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
  timer.stop();
  std::cerr << "Remeshing: " << timer.time() << '\n';

  timer.reset();
  timer.start();
  // Output the facets of the c3t3 to an OFF file. The facets will not be
  // oriented.
  std::ofstream off_file("out.off");
  c3t3.output_boundary_to_off(off_file);
  timer.stop();
  std::cerr << "Export surface: " << timer.time() << '\n';

  return off_file.fail() ? EXIT_FAILURE : EXIT_SUCCESS;
}

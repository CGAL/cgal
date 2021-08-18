#define CGAL_MESH_3_VERBOSE 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>

// Domain
using K           = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron  = CGAL::Polyhedron_3<K>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_3<Polyhedron, K>;

using Concurrency_tag =
#ifdef CGAL_CONCURRENT_MESH_3
  CGAL::Parallel_tag;
#else
  CGAL::Sequential_tag;
#endif

// Triangulation
using Tr   = CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

// Criteria
using Facet_SI_visitor = CGAL::Mesh_3::Facet_criterion_visitor_with_self_intersections<Tr>;
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr,
                              CGAL::Mesh_edge_criteria_3<Tr>,
                              CGAL::Mesh_facet_criteria_3<Tr, Facet_SI_visitor>,
                              CGAL::Mesh_cell_criteria_3<Tr> >;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char*argv[])
{
  const char* fname = (argc>1)?argv[1]:"data/elephant.off";
  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input(fname);
  input >> polyhedron;
  if(input.fail()){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  input.close();

  if (!CGAL::is_triangle_mesh(polyhedron)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // Self-intersections
  typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tri_pairs;
  PMP::self_intersections<Concurrency_tag>(
    faces(polyhedron),
    polyhedron,
    std::back_inserter(intersected_tri_pairs));
  std::cout << intersected_tri_pairs.size() << " pairs of triangles intersect." << std::endl;

  using FProperty_tag = CGAL::dynamic_face_property_t<bool>;
  using SIMap         = boost::property_map<Polyhedron, FProperty_tag>::type;
  SIMap sif_pmap = get(CGAL::dynamic_face_property_t<bool>(), polyhedron);

  for (face_descriptor f : faces(polyhedron))
    put(sif_pmap, f, false);

  for (std::pair<face_descriptor, face_descriptor> p : intersected_tri_pairs)
  {
    put(sif_pmap, p.first, true);
    put(sif_pmap, p.second, true);
  }

  // Create domain
  Mesh_domain domain(polyhedron, nullptr, sif_pmap);

  CGAL::Mesh_facet_criteria_3<Tr, Facet_SI_visitor> facet_criteria(
    25.,//angle_bound
    0.15,//radius bound
    0.003);//distance bound
  CGAL::Mesh_cell_criteria_3<Tr> cell_criteria(
    3.,
    0.);//radius_edge_bound

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  // Output
  CGAL::dump_c3t3(c3t3, "out_self_intersections ");
//  std::ofstream medit_file("out_1.mesh");
//  c3t3.output_to_medit(medit_file);
//  medit_file.close();
//
//  // Set tetrahedron size (keep cell_radius_edge_ratio), ignore facets
//  Mesh_criteria new_criteria(cell_radius_edge_ratio=3, cell_size=0.03);
//
//  // Mesh refinement (and make the output manifold)
//  CGAL::refine_mesh_3(c3t3, domain, new_criteria, manifold());
//
//  // Output
//  medit_file.open("out_2.mesh");
//  c3t3.output_to_medit(medit_file);

  return EXIT_SUCCESS;
}

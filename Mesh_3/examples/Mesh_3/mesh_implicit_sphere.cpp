#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Triangulation_lazy_ds_cell_base_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Implicit_mesh_domain_3<Function,K> Mesh_domain;

// Triangulation
#if defined(CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE) \
 || defined(CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE)
  typedef CGAL::Kernel_traits<Mesh_domain>::Kernel                        PMDKernel;
  typedef CGAL::details::Mesh_geom_traits_generator<PMDKernel>::type      Geom_traits;
  typedef CGAL::Triangulation_lazy_ds_cell_base_3<>                       DS_cell_base;
  typedef CGAL::Triangulation_cell_base_with_circumcenter_3<
            Geom_traits, DS_cell_base>                                    Cell_base_with_cc;
  typedef CGAL::Regular_triangulation_cell_base_3<
            Geom_traits, Cell_base_with_cc>                               Regular_cell_base;
  typedef CGAL::Mesh_triangulation_3<
              Mesh_domain,
              K,
              Geom_traits,
              Regular_cell_base>::type                                    Tr;
#else
  typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Function
FT sphere_function (const Point& p)
{ return CGAL::squared_distance(p, Point(CGAL::ORIGIN))-1; }

int main()
{
  // Domain (Warning: Sphere_3 constructor uses squared radius !)
  Mesh_domain domain(sphere_function, K::Sphere_3(CGAL::ORIGIN, 2.));

  // Mesh criteria
  Mesh_criteria criteria(facet_angle=30, facet_size=0.1, facet_distance=0.025,
                         cell_radius_edge_ratio=2, cell_size=0.1);
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);

  return 0;
}


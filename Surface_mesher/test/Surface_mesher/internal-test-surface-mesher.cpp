#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits Gt;
typedef Gt::Sphere_3 Sphere_3;

typedef double (*Function)(double, double, double);
typedef CGAL::Implicit_surface_3<Gt, Function> Surface_3;
typedef CGAL::Surface_mesh_traits_generator_3<Surface_3>::type SMTraits;

typedef CGAL::Surface_mesh_default_criteria_3<Tr> Criteria;

// shorter names
using CGAL::Surface_mesher::Surface_mesher;
using CGAL::Surface_mesher::Surface_mesher_base;
using CGAL::Surface_mesher::Surface_mesher_manifold_base;
using CGAL::Surface_mesher::Surface_mesher_regular_edges_base;
using CGAL::Surface_mesher::Surface_mesher_regular_edges_without_boundary_base;

// define the surface meshers
typedef Surface_mesher_base<C2t3, Surface_3, SMTraits, Criteria> SM_base;
typedef Surface_mesher_regular_edges_base<C2t3, Surface_3, SMTraits, Criteria> SMRE_base;
typedef Surface_mesher_regular_edges_without_boundary_base<C2t3, Surface_3, SMTraits, Criteria> SMREWB_base;
typedef Surface_mesher_manifold_base<C2t3, Surface_3, SMTraits, Criteria, SMRE_base> SMM_base;
typedef Surface_mesher_manifold_base<C2t3, Surface_3, SMTraits, Criteria, SMREWB_base> SMMWB_base;

typedef Surface_mesher<SM_base> SM;
typedef Surface_mesher<SMRE_base> SMRE;
typedef Surface_mesher<SMREWB_base> SMREWB;
typedef Surface_mesher<SMM_base> SMM;
typedef Surface_mesher<SMMWB_base> SMWB;

// typedef SM Mesher;
// typedef SMRE Mesher;
typedef SMREWB Mesher;
// typedef SMM Mesher;
// typedef SMMWB Mesher;

double sphere_function (double x, double y, double z) {
  const double x2=x*x, y2=y*y, z2=z*z;
  return x2+y2+z2-1;
}

int main(int, char **) {
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // defining the surface
  Surface_3 surface(sphere_function,            // pointer to function
                    Sphere_3(CGAL::ORIGIN, 2.), // bounding sphere
                    1e-03);  // precision for intersections computations

  // defining meshing criteria
  Criteria criteria(30.,  //angular bound
                    0.1,  //radius bound
                    0.1); //distance bound

  SMTraits().construct_initial_points_object()(surface,
                                               CGAL::inserter(tr),
                                               20);

  Mesher mesher(c2t3, surface, SMTraits(), criteria);
  mesher.refine_mesh(true); // true == verbosity

  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}

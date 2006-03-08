#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/MP_Float.h>
#include <CGAL/CORE_Expr.h>

#include <CGAL/make_surface_mesh.h>

#define CGAL_SURFACE_MESHER_TEST 1
#include <CGAL/Surface_mesh_triangulation_generator_3.h> // undocumented

#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Implicit_surface_3.h>

#include <iostream>

struct Sphere {
  template <typename FT>
  FT operator()(FT x, FT y, FT z)
  {
    return x*x+y*y+z*z-1.;
  }
};

enum Flag { NO_NOTHING, DO_NOT_RUN };

template <typename K>
void test_with_kernel(K, Flag flag = NO_NOTHING)
{
  typedef typename CGAL::Surface_mesh_triangulation_generator_3<K>::Type Tr;
  typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2T3;

  typedef typename K::Sphere_3 Sphere_3;

  using CGAL::make_surface_mesh;

  Tr tr;

  double angle_bound = 30.;
  double radius_bound = 0.1;
  double distance_bound = 0.1;
  int initial_number_of_points = 20;

  if(flag == DO_NOT_RUN)
  {
    angle_bound = radius_bound = distance_bound = 1e6;
    initial_number_of_points = 0;
  }

  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle_bound,
                                                     radius_bound,
                                                     distance_bound);

  // 2D-complex in 3D-Delaunay triangulation
  C2T3 c2t3 (tr);

  // Surface meshing
  make_surface_mesh(c2t3,
		    CGAL::make_implicit_surface_3(K(),
						  Sphere(),
						  Sphere_3(CGAL::ORIGIN, 2.),
						  1e-03),
		    criteria,
		    CGAL::Non_manifold_tag(),
                    initial_number_of_points);
  
  std::cout << "Final number of points: " << tr.number_of_vertices() 
            << std::endl;
}

int main(int argc, char **argv)
{
  test_with_kernel(CGAL::Exact_predicates_inexact_constructions_kernel());

  test_with_kernel(CGAL::Filtered_kernel<CGAL::Cartesian<float> >());

  test_with_kernel(CGAL::Filtered_kernel<
                   CGAL::Simple_cartesian<CORE::Expr> >(),
                   DO_NOT_RUN);

  test_with_kernel(CGAL::Cartesian<CGAL::Lazy_exact_nt<double> >(),
                   DO_NOT_RUN);
}

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

#include <CGAL/Timer.h>

#include <iostream>

struct Sphere {
  template <typename FT>
  FT operator()(FT x, FT y, FT z)
  {
    return x*x+y*y+z*z-FT(1);
  }
};

enum Flag { DEFAULT, DO_NOT_RUN};

template <typename K, typename Tag>
struct Test_with_kernel {
  void operator()(Flag flag = DEFAULT)
  {
    if(flag == DO_NOT_RUN) return;

    typedef typename CGAL::Surface_mesh_triangulation_generator_3<K>::Type Tr;
    typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2T3;

    typedef typename K::Sphere_3 Sphere_3;

    using CGAL::make_surface_mesh;

    Tr tr;

    const double angle_bound = 30.;
    const double radius_bound = 0.1;
    const double distance_bound = 0.1;
    const int initial_number_of_points = 20;

    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle_bound,
                                                       radius_bound,
                                                       distance_bound);

    std::cout << "   implicit_surface sphere(O, 1.)\n";

    // 2D-complex in 3D-Delaunay triangulation
    C2T3 c2t3 (tr);

    CGAL::Timer timer;

    timer.start();
    // Surface meshing
    make_surface_mesh(c2t3,
                      CGAL::make_implicit_surface_3(K(),
                                                    Sphere(),
                                                    Sphere_3(CGAL::ORIGIN, 4.),
                                                    1e-03),
                      criteria,
                      Tag(),
                      initial_number_of_points);
    timer.stop();

    std::cout << "Final number of points: " << tr.number_of_vertices()
              << "  (elasped time: " << timer.time() << ")\n\n";

    // same test, with a Sphere_3
    std::cout << "   Kernel::Sphere_3(ORIGIN, 1.)\n";
    Tr tr_2;
    C2T3 c2t3_2(tr_2);
    timer.reset(); timer.start();
    make_surface_mesh(c2t3_2,
                      Sphere_3(CGAL::ORIGIN, 1.),
                      criteria,
                      Tag(),
                      initial_number_of_points);  
    timer.stop();
    std::cout << "Final number of points: " << tr_2.number_of_vertices() 
              << "  (elasped time: " << timer.time() << ")\n\n";
  
  }
};

template <typename Tag>
void test_with_tag(Tag = CGAL::Non_manifold_tag())
{
  std::cout << "\nKERNEL "
    "CGAL::Exact_predicates_inexact_constructions_kernel...\n";

  Test_with_kernel<CGAL::Exact_predicates_inexact_constructions_kernel,
    Tag>()();

  std::cout << "\nKERNEL CGAL::Filtered_kernel<CGAL::Cartesian<float> >...\n";
  Test_with_kernel<CGAL::Filtered_kernel<CGAL::Cartesian<float> >,Tag>()();

  Test_with_kernel<CGAL::Filtered_kernel<CGAL::Simple_cartesian<CORE::Expr> >,
    Tag>()(DO_NOT_RUN);

  Test_with_kernel<CGAL::Cartesian<CGAL::Lazy_exact_nt<double> >,
    Tag >()(DO_NOT_RUN);
}

int main(int, char **)
{
  std::cout << "\n\n    NON MANIFOLD VERSION...\n";
  test_with_tag(CGAL::Non_manifold_tag());
  std::cout << "\n\n    MANIFOLD WITH BOUNDARY VERSION...\n";
  test_with_tag(CGAL::Manifold_with_boundary_tag());
  std::cout << "\n\n    MANIFOLD VERSION...\n";
  test_with_tag(CGAL::Manifold_tag());
}

// TODO: use this C++ trick:
// // Explicit instantiation of the whole class :
// template class CGAL::Triangulation_2<TestK>;


#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/make_surface_mesh.h>

#define CGAL_SURFACE_MESHER_TEST 1
#include <CGAL/Surface_mesh_triangulation_generator_3.h> // undocumented

#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/Timer.h>

#include <iostream>

template <typename K>
class Sphere {
public:
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Sphere_3 Sphere_3;

  Sphere(Sphere_3 sphere)
    : sphere(sphere)
  {
  }

  FT operator()(const Point_3& p) const
  {
    FT x = p.x();
    FT y = p.y();
    FT z = p.z();
    x-=sphere.center().x();
    y-=sphere.center().y();
    z-=sphere.center().z();
    return x*x+y*y+z*z-sphere.squared_radius();
  }
private:
  Sphere_3 sphere;
};

template <typename K>
class Two_spheres
{
public:
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;

  Two_spheres(Sphere<K> sphere1, Sphere<K> sphere2)
    :  sphere1(sphere1), sphere2(sphere2)
  {
  }

  FT operator()(const Point_3 p) const
  {
    return sphere1(p)*sphere2(p);
  }
private:
  Sphere<K> sphere1;
  Sphere<K> sphere2;
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
    typedef typename K::Point_3 Point_3;

    using CGAL::make_surface_mesh;

    Tr tr;

    const double angle_bound = 30.;
    const double radius_bound = 0.1;
    const double distance_bound = 0.1;
    const int initial_number_of_points = 20;

    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle_bound,
                                                       radius_bound,
                                                       distance_bound);

    std::cout << "   implicit_surface: sphere\n";

    // 2D-complex in 3D-Delaunay triangulation
    C2T3 c2t3 (tr);

    typedef CGAL::Implicit_surface_3<K, Sphere<K> > Surface;

    CGAL::Timer timer;

    timer.start();
    // Surface meshing
    CGAL::make_surface_mesh(c2t3,
                            Surface(Sphere<K>(Sphere_3(Point_3(0.3, -5., 1/3.),
                                                       1.)),
                                    Sphere_3(Point_3(0.1, -4.5, 0.), 3.*3.),
                                    1e-03),
                            criteria,
                            Tag(),
                            initial_number_of_points);
    timer.stop();

    std::cout << "Final number of points: " << tr.number_of_vertices()
              << "  (elapsed time: " << timer.time() << ")\n\n";

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
              << "  (elapsed time: " << timer.time() << ")\n\n";

    typedef CGAL::Implicit_surface_3<K, Two_spheres<K> > Surface2;
    typedef typename CGAL::Surface_mesh_traits_generator_3<Surface>::Type Surface_mesh_traits;
    // initial points for a Sphere<K>

    // same test, with Two_spheres
    std::cout << "   Two spheres\n";
    Tr tr_3;

    Sphere<K> sphere1 = Sphere<K>(Sphere_3(CGAL::ORIGIN, 1.));
    Sphere<K> sphere2 = Sphere<K>(Sphere_3(Point_3(0.5, 0., 0.), 0.49*0.49));

    // trick: insert in tr_3 the initial points for sphere2
    Surface surface_of_sphere_2(sphere2, Sphere_3(Point_3(0.5, 0., 0.), 2.));

    Surface_mesh_traits().construct_initial_points_object()
      (surface_of_sphere_2, CGAL::inserter(tr_3), initial_number_of_points);

    C2T3 c2t3_3(tr_3);
    timer.reset(); timer.start();
    make_surface_mesh(c2t3_3,
                      Surface2(Two_spheres<K>(sphere1,sphere2),
                               Sphere_3(CGAL::ORIGIN, 2.)),
                      criteria,
                      Tag(),
                      initial_number_of_points);
    timer.stop();
    std::cout << "Final number of points: " << tr_3.number_of_vertices()
              << "  (elapsed time: " << timer.time() << ")\n\n";


  }
};

template <typename Tag>
void test_with_tag(Tag = CGAL::Non_manifold_tag())
{
  std::cout << "\nKERNEL "
    "CGAL::Exact_predicates_inexact_constructions_kernel...\n";

  Test_with_kernel<CGAL::Exact_predicates_inexact_constructions_kernel,
    Tag>()();

#ifndef CGAL_SURFACE_MESHER_SINGLE_TEST
  std::cout << "\nKERNEL CGAL::Filtered_kernel<CGAL::Cartesian<double> >...\n";
  Test_with_kernel<CGAL::Filtered_kernel<CGAL::Cartesian<double> >,Tag>()();

  Test_with_kernel<CGAL::Cartesian<CGAL::Lazy_exact_nt<double> >,
    Tag >()(DO_NOT_RUN);
#endif // #ifndef CGAL_SURFACE_MESHER_SINGLE_TEST
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


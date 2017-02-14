#define CGAL_SURFACE_MESHER_VERBOSE
//#define CGAL_SURFACE_MESHER_DEBUG_INITIAL_POINTS
//#define CGAL_MESHES_DEBUG_REFINEMENT_POINTS
#include <CGAL/basic.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/MP_Float.h>
#ifdef CGAL_USE_CORE
#  include <CGAL/CORE_Expr.h>
#endif

#include <CGAL/Surface_mesh_periodic_triangulation_3.h> 
#include <CGAL/Surface_mesh_periodic_criteria_3.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/make_surface_mesh.h>

#include <CGAL/Timer.h>

#include <iostream>

template <typename K>
class Schwarz_p {
public:
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;

  Schwarz_p(Vector_3 v, FT diff = FT(1) )
    : v(v), diff(diff)
  {
  }

  FT operator()(const Point_3& p) const
  {
    Point_3 p2 = p+v;
    FT x = std::cos(diff*p2.x()*2*std::acos(-1));
    FT y = std::cos(diff*p2.y()*2*std::acos(-1));
    FT z = std::cos(diff*p2.z()*2*std::acos(-1));
    return x+y+z;
  }
private:
  Vector_3 v;
  FT diff;
};

template <typename K>
class Double_p 
{
public:
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  
  Double_p(Vector_3 v, FT diff)
    :  sp1(v), sp2(v,diff)
  {
  }

  FT operator()(const Point_3 p) const
  {
    return sp1(p)*sp2(p);
  }
private:
  Schwarz_p<K> sp1;
  Schwarz_p<K> sp2;
};
  

enum Flag { DEFAULT, DO_NOT_RUN};

template <typename K, typename Tag>
struct Test_with_kernel {
  void operator()(Flag flag = DEFAULT)
  {
    if(flag == DO_NOT_RUN) return;

    typedef CGAL::Surface_mesh_periodic_triangulation_3 Tr;
    typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2T3;

    typedef typename K::Point_3 Point_3;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Sphere_3 Sphere_3;
    
    using CGAL::make_surface_mesh;

    Tr tr;

    const double angle_bound = 30.;
    const double radius_bound = 0.1;
    const double distance_bound = 0.1;
    const int initial_number_of_points = 20;

    CGAL::Surface_mesh_periodic_criteria_3<Tr> criteria(&tr, 
	angle_bound, radius_bound, distance_bound);

    std::cout << "   implicit_surface: schwarz_p\n";

    // 2D-complex in 3D-Delaunay triangulation
    C2T3 c2t3 (tr);

    typedef CGAL::Implicit_surface_3<K, Schwarz_p<K> > Surface;

    CGAL::Timer timer;

    timer.start();
    // Surface meshing
    Surface surface(Schwarz_p<K>(Vector_3(0.3, -5., 1/3.),1.),
	Sphere_3(Point_3(0.1, -4.5, 0.), 2.),
	1e-03);

    CGAL::make_surface_mesh(c2t3,
	surface,
	criteria,
	Tag(),
	initial_number_of_points);
    timer.stop();

    std::cout << "Final number of points: " << tr.number_of_vertices()
              << "  (elasped time: " << timer.time() << ")\n\n";

    typedef CGAL::Implicit_surface_3<K, Double_p<K> > Surface2;
    typedef typename CGAL::Surface_mesh_traits_generator_3<Surface>::Type
      Surface_mesh_traits;
    // initial points for Double_p<K>

    // same test, with Two_spheres
    std::cout << "   Double_p\n";
    Tr tr_2;

    Schwarz_p<K> sp1 = Schwarz_p<K>(Vector_3(0.1,-4.5,0.));

    // trick: insert in tr_2 the initial points for sp2
    Surface surface_of_sp1(sp1, Sphere_3(CGAL::ORIGIN, 2.));

    Surface_mesh_traits().construct_initial_points_object()
      (surface_of_sp1, CGAL::inserter(tr_2), initial_number_of_points);

    C2T3 c2t3_2(tr_2);
    timer.reset(); timer.start();
    make_surface_mesh(c2t3_2,
	Surface2(Double_p<K>(Vector_3(0.1,-4.5,0.),2),
	    Sphere_3(CGAL::ORIGIN, 2.)),
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


  //TODO
// #ifndef CGAL_SURFACE_MESHER_SINGLE_TEST
//   std::cout << "\nKERNEL CGAL::Filtered_kernel<CGAL::Cartesian<float> >...\n";
//   Test_with_kernel<CGAL::Filtered_kernel<CGAL::Cartesian<float> >,Tag>()();

//   Test_with_kernel<CGAL::Cartesian<CGAL::Lazy_exact_nt<double> >,
//     Tag >()(DO_NOT_RUN);
// #endif // #ifndef CGAL_SURFACE_MESHER_SINGLE_TEST
}

int main(int, char **)
{
  std::cout << "\n\n    NON MANIFOLD VERSION...\n";
  test_with_tag(CGAL::Non_manifold_tag());
  //TODO
//   std::cout << "\n\n    MANIFOLD WITH BOUNDARY VERSION...\n";
//   test_with_tag(CGAL::Manifold_with_boundary_tag());
//   std::cout << "\n\n    MANIFOLD VERSION...\n";
//   test_with_tag(CGAL::Manifold_tag());
}

// TODO: use this C++ trick:
// // Explicit instantiation of the whole class :
// template class CGAL::Triangulation_2<TestK>;


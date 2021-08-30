#include <CGAL/use.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// This just tests if every call we promise is actually available
typedef CGAL::Exact_predicates_exact_constructions_kernel Epic;

typedef CGAL::Cartesian<double>     Cartesian;
typedef Epic K;
typedef CGAL::Point_3< K >          P;
typedef CGAL::Segment_3< K >        S;
typedef CGAL::Line_3< K >           L;
typedef CGAL::Plane_3< K >          Pl;
typedef CGAL::Triangle_3< K >       Tr;
typedef CGAL::Ray_3< K >            R;
typedef CGAL::Iso_cuboid_3< K >     Cub;

typedef CGAL::Sphere_3< K >         Sph;
typedef CGAL::Circle_3< K >         C;

typedef CGAL::Tetrahedron_3< K >    T;

typedef CGAL::Bbox_3                Bbox_3;

template<class A, class B>
void call_intersection_global(const A& a, const B& b) {
  const auto x = CGAL::intersection(a, b);
  const auto y = CGAL::intersection(b, a);
  const auto z = CGAL::intersection(b, a);
  CGAL_USE(x);
  CGAL_USE(y);
  CGAL_USE(z);
}

template<class A, class B>
void call_do_intersect_global(const A& a, const B& b) {
  CGAL::do_intersect(a, b);
  CGAL::do_intersect(b, a);
}

template<class A, class B, class K>
void call_intersection_with_kernel(const A& a, const B& b, const K&) {
  typedef typename K::Intersect_3 Intersect;
  const auto x = Intersect()(a, b);
  const auto y = Intersect()(b, a);
}

template<class A, class B, class K>
void call_do_intersect_with_kernel(const A& a, const B& b, const K&) {
  typedef typename K::Do_intersect_3 Do_inter;
  Do_inter()(a, b);
  Do_inter()(b, a);
}


int main(int argc, char**)
{
  CGAL::Interval_nt_advanced::Protector p;
  CGAL_USE(p);
  //we only want to check compilation
  if(argc > 666 )
  {
    call_intersection_global(S(), S());
    call_intersection_global(S(), L());
    call_intersection_global(S(), Pl());
    call_intersection_global(S(), Tr());
    call_intersection_global(S(), R());
    call_intersection_global(S(), Cub());

    call_intersection_global(L(), S());
    call_intersection_global(L(), L());
    call_intersection_global(L(), Pl());
    call_intersection_global(L(), Tr());
    call_intersection_global(L(), R());
    call_intersection_global(L(), Cub());

    call_intersection_global(Pl(), S());
    call_intersection_global(Pl(), L());
    call_intersection_global(Pl(), Pl());
    call_intersection_global(Pl(), Tr());
    call_intersection_global(Pl(), R());
    // call_intersection_global(Pl(), Cub());

    // special
    const auto plplpl = CGAL::intersection(Pl(), Pl(), Pl());

    call_intersection_global(Tr(), S());
    call_intersection_global(Tr(), L());
    call_intersection_global(Tr(), Pl());
    call_intersection_global(Tr(), Tr());
    call_intersection_global(Tr(), R());
    // call_intersection_global(Tr(), Cub());

    call_intersection_global(R(), S());
    call_intersection_global(R(), L());
    call_intersection_global(R(), Pl());
    call_intersection_global(R(), Tr());
    call_intersection_global(R(), R());
    call_intersection_global(R(), Cub());

    call_intersection_global(Cub(), S());
    call_intersection_global(Cub(), L());
    // call_intersection_global(Cub(), Pl());
    // call_intersection_global(Cub(), Tr());
    call_intersection_global(Cub(), R());
    call_intersection_global(Cub(), Cub());
//    call_intersection_global(Cub(), Bbox_3());

    call_intersection_global(Bbox_3(), L());
    call_intersection_global(Bbox_3(), S());
    call_intersection_global(Bbox_3(), R());
    call_intersection_global(Bbox_3(), Cub());
    CGAL::intersection(Bbox_3(), Bbox_3());

    call_intersection_global(T(), L());


    // with kernel

    call_intersection_with_kernel(S(), S(), K());
    call_intersection_with_kernel(S(), L(), K());
    call_intersection_with_kernel(S(), Pl(), K());
    call_intersection_with_kernel(S(), Tr(), K());
    call_intersection_with_kernel(S(), R(), K());
    call_intersection_with_kernel(S(), Cub(), K());

    call_intersection_with_kernel(L(), S(), K());
    call_intersection_with_kernel(L(), L(), K());
    call_intersection_with_kernel(L(), Pl(), K());
    call_intersection_with_kernel(L(), Tr(), K());
    call_intersection_with_kernel(L(), R(), K());
    call_intersection_with_kernel(L(), Cub(), K());

    call_intersection_with_kernel(Pl(), S(), K());
    call_intersection_with_kernel(Pl(), L(), K());
    call_intersection_with_kernel(Pl(), Pl(), K());
    call_intersection_with_kernel(Pl(), Tr(), K());
    call_intersection_with_kernel(Pl(), R(), K());
    // call_intersection_with_kernel(Pl(), Cub(), K());

    //special
    K::Intersect_3()(Pl(), Pl(), Pl());

    call_intersection_with_kernel(Tr(), S(), K());
    call_intersection_with_kernel(Tr(), L(), K());
    call_intersection_with_kernel(Tr(), Pl(), K());
    call_intersection_with_kernel(Tr(), Tr(), K());
    call_intersection_with_kernel(Tr(), R(), K());
    // call_intersection_with_kernel(Tr(), Cub(), K());

    call_intersection_with_kernel(R(), S(), K());
    call_intersection_with_kernel(R(), L(), K());
    call_intersection_with_kernel(R(), Pl(), K());
    call_intersection_with_kernel(R(), Tr(), K());
    call_intersection_with_kernel(R(), R(), K());
    call_intersection_with_kernel(R(), Cub(), K());

    call_intersection_with_kernel(Cub(), S(), K());
    call_intersection_with_kernel(Cub(), L(), K());
    // call_intersection_with_kernel(Cub(), Pl(), K());
    // call_intersection_with_kernel(Cub(), Tr(), K());
    call_intersection_with_kernel(Cub(), R(), K());
    call_intersection_with_kernel(Cub(), Cub(), K());

    call_intersection_with_kernel(Bbox_3(), L(), K());
    call_intersection_with_kernel(Bbox_3(), S(), K());
    call_intersection_with_kernel(Bbox_3(), R(), K());

    // The doc defines calls to do_intersect for these objects

    // Plane_3<Kernel>
    // Line_3<Kernel>
    // Ray_3<Kernel>
    // Segment_3<Kernel>
    // Triangle_3<Kernel>.
    // Bbox_3.
    call_do_intersect_global(Pl(), Pl());
    call_do_intersect_global(Pl(), L());
    call_do_intersect_global(Pl(), R());
    call_do_intersect_global(Pl(), S());
    call_do_intersect_global(Pl(), Tr());
    call_do_intersect_global(Pl(), Bbox_3());

    call_do_intersect_global(L(), Pl());
    call_do_intersect_global(L(), L());
    call_do_intersect_global(L(), R());
    call_do_intersect_global(L(), S());
    call_do_intersect_global(L(), Tr());
    call_do_intersect_global(L(), Bbox_3());

    call_do_intersect_global(R(), Pl());
    call_do_intersect_global(R(), L());
    call_do_intersect_global(R(), R());
    call_do_intersect_global(R(), S());
    call_do_intersect_global(R(), Tr());
    call_do_intersect_global(R(), Bbox_3());

    call_do_intersect_global(S(), Pl());
    call_do_intersect_global(S(), L());
    call_do_intersect_global(S(), R());
    call_do_intersect_global(S(), S());
    call_do_intersect_global(S(), Tr());
    call_do_intersect_global(S(), Bbox_3());

    call_do_intersect_global(Tr(), Pl());
    call_do_intersect_global(Tr(), L());
    call_do_intersect_global(Tr(), R());
    call_do_intersect_global(Tr(), S());
    call_do_intersect_global(Tr(), Tr());
    call_do_intersect_global(Tr(), Bbox_3());

    call_do_intersect_global(Tr(), Pl());
    call_do_intersect_global(Tr(), L());
    call_do_intersect_global(Tr(), R());
    call_do_intersect_global(Tr(), S());
    call_do_intersect_global(Tr(), Tr());
    call_do_intersect_global(Tr(), Bbox_3());

    call_do_intersect_global(Bbox_3(), Pl());
    call_do_intersect_global(Bbox_3(), L());
    call_do_intersect_global(Bbox_3(), R());
    call_do_intersect_global(Bbox_3(), S());
    call_do_intersect_global(Bbox_3(), Tr());
    call_do_intersect_global(Bbox_3(), Sph());
    call_do_intersect_global(Bbox_3(), Bbox_3());

    // with_kernel
    call_do_intersect_with_kernel(Pl(), Pl(), K());
    call_do_intersect_with_kernel(Pl(), L(), K());
    call_do_intersect_with_kernel(Pl(), R(), K());
    call_do_intersect_with_kernel(Pl(), S(), K());
    call_do_intersect_with_kernel(Pl(), Tr(), K());
    call_do_intersect_with_kernel(Pl(), Bbox_3(), K());

    call_do_intersect_with_kernel(L(), Pl(), K());
    call_do_intersect_with_kernel(L(), L(), K());
    call_do_intersect_with_kernel(L(), R(), K());
    call_do_intersect_with_kernel(L(), S(), K());
    call_do_intersect_with_kernel(L(), Tr(), K());
    call_do_intersect_with_kernel(L(), Bbox_3(), K());

    call_do_intersect_with_kernel(R(), Pl(), K());
    call_do_intersect_with_kernel(R(), L(), K());
    call_do_intersect_with_kernel(R(), R(), K());
    call_do_intersect_with_kernel(R(), S(), K());
    call_do_intersect_with_kernel(R(), Tr(), K());
    call_do_intersect_with_kernel(R(), Bbox_3(), K());

    call_do_intersect_with_kernel(S(), Pl(), K());
    call_do_intersect_with_kernel(S(), L(), K());
    call_do_intersect_with_kernel(S(), R(), K());
    call_do_intersect_with_kernel(S(), S(), K());
    call_do_intersect_with_kernel(S(), Tr(), K());
    call_do_intersect_with_kernel(S(), Bbox_3(), K());

    call_do_intersect_with_kernel(Tr(), Pl(), K());
    call_do_intersect_with_kernel(Tr(), L(), K());
    call_do_intersect_with_kernel(Tr(), R(), K());
    call_do_intersect_with_kernel(Tr(), S(), K());
    call_do_intersect_with_kernel(Tr(), Tr(), K());
    call_do_intersect_with_kernel(Tr(), Bbox_3(), K());

    call_do_intersect_with_kernel(Tr(), Pl(), K());
    call_do_intersect_with_kernel(Tr(), L(), K());
    call_do_intersect_with_kernel(Tr(), R(), K());
    call_do_intersect_with_kernel(Tr(), S(), K());
    call_do_intersect_with_kernel(Tr(), Tr(), K());
    call_do_intersect_with_kernel(Tr(), Bbox_3(), K());

    call_do_intersect_with_kernel(Bbox_3(), Pl(), K());
    call_do_intersect_with_kernel(Bbox_3(), L(), K());
    call_do_intersect_with_kernel(Bbox_3(), R(), K());
    call_do_intersect_with_kernel(Bbox_3(), S(), K());
    call_do_intersect_with_kernel(Bbox_3(), Sph(), K());
    call_do_intersect_with_kernel(Bbox_3(), Tr(), K());
  }
  return EXIT_SUCCESS;
}

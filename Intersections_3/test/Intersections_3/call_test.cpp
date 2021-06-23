#include <CGAL/use.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/intersections.h>

// This just tests if every call we promise is actually available

template<class A, class B>
void call_intersection_global(const A& a, const B& b)
{
  const auto x = CGAL::intersection(a, b);
  const auto y = CGAL::intersection(b, a);
  const auto z = CGAL::intersection(b, a);
  CGAL_USE(x);
  CGAL_USE(y);
  CGAL_USE(z);
}

template<class A, class B>
void call_do_intersect_global(const A& a, const B& b)
{
  const auto x = CGAL::do_intersect(a, b);
  const auto y = CGAL::do_intersect(b, a);
  CGAL_USE(x);
  CGAL_USE(y);
}

template<class A, class B, class K>
void call_intersection_with_kernel(const A& a, const B& b, const K& k)
{
  typename K::Intersect_3 intersect = k.intersect_3_object();
  const auto x = intersect(a, b);
  const auto y = intersect(b, a);
  CGAL_USE(x);
  CGAL_USE(y);
}

template<class A, class B, class K>
void call_do_intersect_with_kernel(const A& a, const B& b, const K& k)
{
  typename K::Do_intersect_3 do_intersect = k.do_intersect_3_object();
  const bool x = do_intersect(a, b);
  const bool y = do_intersect(b, a);
  CGAL_USE(x);
  CGAL_USE(y);
}

template <typename K>
void test(const int argc)
{
  typedef CGAL::Iso_cuboid_3< K >     Cub;
  typedef CGAL::Line_3< K >           L;
  typedef CGAL::Plane_3< K >          Pl;
  typedef CGAL::Point_3< K >          P;
  typedef CGAL::Ray_3< K >            R;
  typedef CGAL::Segment_3< K >        S;
  typedef CGAL::Sphere_3< K >         Sph;
  typedef CGAL::Triangle_3< K >       Tr;
  typedef CGAL::Tetrahedron_3< K >    T;

  typedef CGAL::Bbox_3                Bbox_3;

  //we only want to check compilation
  if(argc > 0)
  {
    // ---------------------------------------------------------------------------------------------
    //                                        INTERSECTION
    // ---------------------------------------------------------------------------------------------

    call_intersection_global(Cub(), Bbox_3());
    call_intersection_global(Cub(), Cub());
    call_intersection_global(Cub(), L());
    call_intersection_global(Cub(), Pl());
    call_intersection_global(Cub(), P());
    call_intersection_global(Cub(), R());
    call_intersection_global(Cub(), S());
    call_intersection_global(Cub(), Tr());

    call_intersection_global(L(), Bbox_3());
    call_intersection_global(L(), Cub());
    call_intersection_global(L(), L());
    call_intersection_global(L(), Pl());
    call_intersection_global(L(), P());
    call_intersection_global(L(), R());
    call_intersection_global(L(), S());
    call_intersection_global(L(), T());
    call_intersection_global(L(), Tr());

    call_intersection_global(Pl(), Bbox_3());
    call_intersection_global(Pl(), Cub());
    call_intersection_global(Pl(), L());
    call_intersection_global(Pl(), Pl());
    call_intersection_global(Pl(), P());
    call_intersection_global(Pl(), R());
    call_intersection_global(Pl(), S());
    call_intersection_global(Pl(), Sph());
    call_intersection_global(Pl(), Tr());
    call_intersection_global(Pl(), T());

    auto plplpl = CGAL::intersection(Pl(), Pl(), Pl());
    CGAL_USE(plplpl);

    call_intersection_global(P(), Bbox_3());
    call_intersection_global(P(), Cub());
    call_intersection_global(P(), L());
    call_intersection_global(P(), Pl());
    call_intersection_global(P(), P());
    call_intersection_global(P(), R());
    call_intersection_global(P(), S());
    call_intersection_global(P(), Sph());
    call_intersection_global(P(), Tr());
    call_intersection_global(P(), T());

    call_intersection_global(R(), Bbox_3());
    call_intersection_global(R(), Cub());
    call_intersection_global(R(), L());
    call_intersection_global(R(), Pl());
    call_intersection_global(R(), P());
    call_intersection_global(R(), R());
    call_intersection_global(R(), S());
    call_intersection_global(R(), Tr());
    call_intersection_global(R(), T());

    call_intersection_global(S(), Bbox_3());
    call_intersection_global(S(), Cub());
    call_intersection_global(S(), L());
    call_intersection_global(S(), Pl());
    call_intersection_global(S(), P());
    call_intersection_global(S(), R());
    call_intersection_global(S(), S());
    call_intersection_global(S(), T());
    call_intersection_global(S(), Tr());

    call_intersection_global(Sph(), Pl());
    call_intersection_global(Sph(), P());
    call_intersection_global(Sph(), Sph());

    call_intersection_global(Tr(), Bbox_3());
    call_intersection_global(Tr(), Cub());
    call_intersection_global(Tr(), L());
    call_intersection_global(Tr(), Pl());
    call_intersection_global(Tr(), P());
    call_intersection_global(Tr(), R());
    call_intersection_global(Tr(), S());
    call_intersection_global(Tr(), Tr());
    call_intersection_global(Tr(), T());

    call_intersection_global(T(), L());
    call_intersection_global(T(), Pl());
    call_intersection_global(T(), P());
    call_intersection_global(T(), R());
    call_intersection_global(T(), S());
    call_intersection_global(T(), Tr());

    const auto bbbb = CGAL::intersection(Bbox_3(), Bbox_3());
    CGAL_USE(bbbb);

    call_intersection_global(Bbox_3(), Cub());
    call_intersection_global(Bbox_3(), L());
    call_intersection_global(Bbox_3(), Pl());
    call_intersection_global(Bbox_3(), P());
    call_intersection_global(Bbox_3(), R());
    call_intersection_global(Bbox_3(), S());
    call_intersection_global(Bbox_3(), Tr());

    call_intersection_global(T(), L());

    // with kernel
    call_intersection_with_kernel(Cub(), Bbox_3(), K());
    call_intersection_with_kernel(Cub(), Cub(), K());
    call_intersection_with_kernel(Cub(), L(), K());
    call_intersection_with_kernel(Cub(), Pl(), K());
    call_intersection_with_kernel(Cub(), P(), K());
    call_intersection_with_kernel(Cub(), R(), K());
    call_intersection_with_kernel(Cub(), S(), K());
    call_intersection_with_kernel(Cub(), Tr(), K());

    call_intersection_with_kernel(L(), Bbox_3(), K());
    call_intersection_with_kernel(L(), Cub(), K());
    call_intersection_with_kernel(L(), L(), K());
    call_intersection_with_kernel(L(), Pl(), K());
    call_intersection_with_kernel(L(), P(), K());
    call_intersection_with_kernel(L(), R(), K());
    call_intersection_with_kernel(L(), S(), K());
    call_intersection_with_kernel(L(), Tr(), K());
    call_intersection_with_kernel(L(), T(), K());

    call_intersection_with_kernel(Pl(), Bbox_3(), K());
    call_intersection_with_kernel(Pl(), Cub(), K());
    call_intersection_with_kernel(Pl(), L(), K());
    call_intersection_with_kernel(Pl(), Pl(), K());
    call_intersection_with_kernel(Pl(), P(), K());
    call_intersection_with_kernel(Pl(), R(), K());
    call_intersection_with_kernel(Pl(), S(), K());
    call_intersection_with_kernel(Pl(), Tr(), K());
    call_intersection_with_kernel(Pl(), T(), K());

    //special
    plplpl = K().intersect_3_object()(Pl(), Pl(), Pl());
    CGAL_USE(plplpl);

    call_intersection_with_kernel(P(), Bbox_3(), K());
    call_intersection_with_kernel(P(), Cub(), K());
    call_intersection_with_kernel(P(), L(), K());
    call_intersection_with_kernel(P(), Pl(), K());
    call_intersection_with_kernel(P(), P(), K());
    call_intersection_with_kernel(P(), R(), K());
    call_intersection_with_kernel(P(), S(), K());
    call_intersection_with_kernel(P(), Sph(), K());
    call_intersection_with_kernel(P(), Tr(), K());
    call_intersection_with_kernel(P(), T(), K());

    call_intersection_with_kernel(R(), Bbox_3(), K());
    call_intersection_with_kernel(R(), Cub(), K());
    call_intersection_with_kernel(R(), L(), K());
    call_intersection_with_kernel(R(), Pl(), K());
    call_intersection_with_kernel(R(), P(), K());
    call_intersection_with_kernel(R(), R(), K());
    call_intersection_with_kernel(R(), S(), K());
    call_intersection_with_kernel(R(), Tr(), K());
    call_intersection_with_kernel(R(), T(), K());

    call_intersection_with_kernel(S(), Bbox_3(), K());
    call_intersection_with_kernel(S(), Cub(), K());
    call_intersection_with_kernel(S(), L(), K());
    call_intersection_with_kernel(S(), Pl(), K());
    call_intersection_with_kernel(S(), P(), K());
    call_intersection_with_kernel(S(), R(), K());
    call_intersection_with_kernel(S(), S(), K());
    call_intersection_with_kernel(S(), Tr(), K());
    call_intersection_with_kernel(S(), T(), K());

    call_intersection_with_kernel(Sph(), Pl(), K());
    call_intersection_with_kernel(Sph(), P(), K());
    call_intersection_with_kernel(Sph(), Sph(), K());

    call_intersection_with_kernel(Tr(), Bbox_3(), K());
    call_intersection_with_kernel(Tr(), Cub(), K());
    call_intersection_with_kernel(Tr(), L(), K());
    call_intersection_with_kernel(Tr(), Pl(), K());
    call_intersection_with_kernel(Tr(), P(), K());
    call_intersection_with_kernel(Tr(), R(), K());
    call_intersection_with_kernel(Tr(), S(), K());
    call_intersection_with_kernel(Tr(), Tr(), K());
    call_intersection_with_kernel(Tr(), T(), K());

    call_intersection_with_kernel(T(), L(), K());
    call_intersection_with_kernel(T(), Pl(), K());
    call_intersection_with_kernel(T(), P(), K());
    call_intersection_with_kernel(T(), R(), K());
    call_intersection_with_kernel(T(), S(), K());
    call_intersection_with_kernel(T(), Tr(), K());

    call_intersection_with_kernel(Bbox_3(), Cub(), K());
    call_intersection_with_kernel(Bbox_3(), L(), K());
    call_intersection_with_kernel(Bbox_3(), Pl(), K());
    call_intersection_with_kernel(Bbox_3(), P(), K());
    call_intersection_with_kernel(Bbox_3(), R(), K());
    call_intersection_with_kernel(Bbox_3(), S(), K());
    call_intersection_with_kernel(Bbox_3(), Tr(), K());

    // ---------------------------------------------------------------------------------------------
    //                                        DO INTERSECT
    // ---------------------------------------------------------------------------------------------

    call_do_intersect_global(L(), Cub());
    call_do_intersect_global(L(), Bbox_3());
    call_do_intersect_global(L(), L());
    call_do_intersect_global(L(), Pl());
    call_do_intersect_global(L(), P());
    call_do_intersect_global(L(), R());
    call_do_intersect_global(L(), S());
    call_do_intersect_global(L(), Sph());
    call_do_intersect_global(L(), Tr());
    call_do_intersect_global(L(), T());

    call_do_intersect_global(Pl(), Bbox_3());
    call_do_intersect_global(Pl(), Cub());
    call_do_intersect_global(Pl(), L());
    call_do_intersect_global(Pl(), Pl());
    call_do_intersect_global(Pl(), P());
    call_do_intersect_global(Pl(), R());
    call_do_intersect_global(Pl(), S());
    call_do_intersect_global(Pl(), Sph());
    call_do_intersect_global(Pl(), Tr());
    call_do_intersect_global(Pl(), T());

    call_do_intersect_global(P(), Bbox_3());
    call_do_intersect_global(P(), Cub());
    call_do_intersect_global(P(), L());
    call_do_intersect_global(P(), Pl());
    call_do_intersect_global(P(), P());
    call_do_intersect_global(P(), R());
    call_do_intersect_global(P(), S());
    call_do_intersect_global(P(), Sph());
    call_do_intersect_global(P(), Tr());
    call_do_intersect_global(P(), T());

    call_do_intersect_global(R(), Bbox_3());
    call_do_intersect_global(R(), Cub());
    call_do_intersect_global(R(), L());
    call_do_intersect_global(R(), Pl());
    call_do_intersect_global(R(), P());
    call_do_intersect_global(R(), R());
    call_do_intersect_global(R(), S());
    call_do_intersect_global(R(), Tr());
    call_do_intersect_global(R(), T());

    call_do_intersect_global(S(), Bbox_3());
    call_do_intersect_global(S(), Cub());
    call_do_intersect_global(S(), Pl());
    call_do_intersect_global(S(), P());
    call_do_intersect_global(S(), L());
    call_do_intersect_global(S(), R());
    call_do_intersect_global(S(), S());
    call_do_intersect_global(S(), Sph());
    call_do_intersect_global(S(), Tr());
    call_do_intersect_global(S(), T());

    call_do_intersect_global(Sph(), Bbox_3());
    call_do_intersect_global(Sph(), Cub());
    call_do_intersect_global(Sph(), Pl());
    call_do_intersect_global(Sph(), P());
    call_do_intersect_global(Sph(), L());
    call_do_intersect_global(Sph(), R());
    call_do_intersect_global(Sph(), S());
    call_do_intersect_global(Sph(), Sph());
    call_do_intersect_global(Sph(), Tr());
    call_do_intersect_global(Sph(), T());

    call_do_intersect_global(Tr(), Bbox_3());
    call_do_intersect_global(Tr(), Cub());
    call_do_intersect_global(Tr(), L());
    call_do_intersect_global(Tr(), Pl());
    call_do_intersect_global(Tr(), P());
    call_do_intersect_global(Tr(), R());
    call_do_intersect_global(Tr(), S());
    call_do_intersect_global(Tr(), Sph());
    call_do_intersect_global(Tr(), Tr());
    call_do_intersect_global(Tr(), T());

    call_do_intersect_global(Tr(), Bbox_3());
    call_do_intersect_global(Tr(), Cub());
    call_do_intersect_global(Tr(), L());
    call_do_intersect_global(Tr(), Pl());
    call_do_intersect_global(Tr(), P());
    call_do_intersect_global(Tr(), R());
    call_do_intersect_global(Tr(), S());
    call_do_intersect_global(Tr(), Sph());
    call_do_intersect_global(Tr(), Tr());
    call_do_intersect_global(Tr(), T());

    call_do_intersect_global(T(), Bbox_3());
    call_do_intersect_global(T(), Cub());
    call_do_intersect_global(T(), L());
    call_do_intersect_global(T(), Pl());
    call_do_intersect_global(T(), P());
    call_do_intersect_global(T(), R());
    call_do_intersect_global(T(), S());
    call_do_intersect_global(T(), Sph());
    call_do_intersect_global(T(), Tr());
    call_do_intersect_global(T(), T());

    call_do_intersect_global(Bbox_3(), Pl());
    call_do_intersect_global(Bbox_3(), L());
    call_do_intersect_global(Bbox_3(), R());
    call_do_intersect_global(Bbox_3(), S());
    call_do_intersect_global(Bbox_3(), Tr());
    call_do_intersect_global(Bbox_3(), Sph());
    call_do_intersect_global(Bbox_3(), Bbox_3());

    // with_kernel
    call_do_intersect_with_kernel(L(), Bbox_3(), K());
    call_do_intersect_with_kernel(L(), Cub(), K());
    call_do_intersect_with_kernel(L(), L(), K());
    call_do_intersect_with_kernel(L(), Pl(), K());
    call_do_intersect_with_kernel(L(), P(), K());
    call_do_intersect_with_kernel(L(), R(), K());
    call_do_intersect_with_kernel(L(), S(), K());
    call_do_intersect_with_kernel(L(), Sph(), K());
    call_do_intersect_with_kernel(L(), Tr(), K());
    call_do_intersect_with_kernel(L(), T(), K());

    call_do_intersect_with_kernel(Pl(), Bbox_3(), K());
    call_do_intersect_with_kernel(Pl(), Cub(), K());
    call_do_intersect_with_kernel(Pl(), L(), K());
    call_do_intersect_with_kernel(Pl(), Pl(), K());
    call_do_intersect_with_kernel(Pl(), P(), K());
    call_do_intersect_with_kernel(Pl(), R(), K());
    call_do_intersect_with_kernel(Pl(), S(), K());
    call_do_intersect_with_kernel(Pl(), Sph(), K());
    call_do_intersect_with_kernel(Pl(), Tr(), K());
    call_do_intersect_with_kernel(Pl(), T(), K());

    call_do_intersect_with_kernel(P(), Bbox_3(), K());
    call_do_intersect_with_kernel(P(), Cub(), K());
    call_do_intersect_with_kernel(P(), L(), K());
    call_do_intersect_with_kernel(P(), Pl(), K());
    call_do_intersect_with_kernel(P(), P(), K());
    call_do_intersect_with_kernel(P(), R(), K());
    call_do_intersect_with_kernel(P(), S(), K());
    call_do_intersect_with_kernel(P(), Sph(), K());
    call_do_intersect_with_kernel(P(), Tr(), K());
    call_do_intersect_with_kernel(P(), T(), K());

    call_do_intersect_with_kernel(R(), Bbox_3(), K());
    call_do_intersect_with_kernel(R(), Cub(), K());
    call_do_intersect_with_kernel(R(), L(), K());
    call_do_intersect_with_kernel(R(), Pl(), K());
    call_do_intersect_with_kernel(R(), P(), K());
    call_do_intersect_with_kernel(R(), R(), K());
    call_do_intersect_with_kernel(R(), S(), K());
    call_do_intersect_with_kernel(R(), Sph(), K());
    call_do_intersect_with_kernel(R(), Tr(), K());
    call_do_intersect_with_kernel(R(), T(), K());

    call_do_intersect_with_kernel(S(), Bbox_3(), K());
    call_do_intersect_with_kernel(S(), Cub(), K());
    call_do_intersect_with_kernel(S(), R(), K());
    call_do_intersect_with_kernel(S(), L(), K());
    call_do_intersect_with_kernel(S(), P(), K());
    call_do_intersect_with_kernel(S(), R(), K());
    call_do_intersect_with_kernel(S(), S(), K());
    call_do_intersect_with_kernel(S(), Sph(), K());
    call_do_intersect_with_kernel(S(), Tr(), K());
    call_do_intersect_with_kernel(S(), T(), K());

    call_do_intersect_with_kernel(Tr(), Bbox_3(), K());
    call_do_intersect_with_kernel(Tr(), Cub(), K());
    call_do_intersect_with_kernel(Tr(), R(), K());
    call_do_intersect_with_kernel(Tr(), L(), K());
    call_do_intersect_with_kernel(Tr(), P(), K());
    call_do_intersect_with_kernel(Tr(), R(), K());
    call_do_intersect_with_kernel(Tr(), S(), K());
    call_do_intersect_with_kernel(Tr(), Tr(), K());
    call_do_intersect_with_kernel(Tr(), T(), K());

    call_do_intersect_with_kernel(T(), Bbox_3(), K());
    call_do_intersect_with_kernel(T(), Cub(), K());
    call_do_intersect_with_kernel(T(), R(), K());
    call_do_intersect_with_kernel(T(), L(), K());
    call_do_intersect_with_kernel(T(), P(), K());
    call_do_intersect_with_kernel(T(), R(), K());
    call_do_intersect_with_kernel(T(), S(), K());
    call_do_intersect_with_kernel(T(), Tr(), K());
    call_do_intersect_with_kernel(T(), T(), K());

    call_do_intersect_with_kernel(Bbox_3(), Cub(), K());
    call_do_intersect_with_kernel(Bbox_3(), L(), K());
    call_do_intersect_with_kernel(Bbox_3(), Pl(), K());
    call_do_intersect_with_kernel(Bbox_3(), P(), K());
    call_do_intersect_with_kernel(Bbox_3(), R(), K());
    call_do_intersect_with_kernel(Bbox_3(), S(), K());
    call_do_intersect_with_kernel(Bbox_3(), Sph(), K());
    call_do_intersect_with_kernel(Bbox_3(), Tr(), K());
  }
}

int main(int argc, char**)
{
  test<CGAL::Cartesian<double> >(argc);
  test<CGAL::Exact_predicates_inexact_constructions_kernel>(argc);
  test<CGAL::Exact_predicates_exact_constructions_kernel>(argc);

  return EXIT_SUCCESS;
}

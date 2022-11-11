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

  // only check compilation
  if(argc > 1)
  {
    // ---------------------------------------------------------------------------------------------
    //                                        INTERSECTION
    // ---------------------------------------------------------------------------------------------

    Bbox_3 bbox { 0, 0, 0, 1, 1, 1 };
    Cub cub { 0, 0, 0, 1, 1, 1 };
    L l { P{0, 0, 0}, P{1, 1, 1} };
    Pl pl { P{0, 0, 0}, P{1, 0, 0}, P{0, 1, 0} };
    P p { 0, 0, 0 };
    R r { P{0, 0, 0}, P{1, 1, 1} };
    S s { P{0, 0, 0}, P{1, 1, 1} };
    Sph sph { P{0, 0, 0}, 1 };
    Tr tr { P{0, 0, 0}, P{1, 0, 0}, P{0, 1, 0} };
    T tet { P{0, 0, 0}, P{1, 0, 0}, P{0, 1, 0}, P{0, 0, 1} };

    call_intersection_global(cub, bbox);
    call_intersection_global(cub, cub);
    call_intersection_global(cub, l);
    call_intersection_global(cub, pl);
    call_intersection_global(cub, p);
    call_intersection_global(cub, r);
    call_intersection_global(cub, s);
    call_intersection_global(cub, tr);

    call_intersection_global(l, bbox);
    call_intersection_global(l, cub);
    call_intersection_global(l, l);
    call_intersection_global(l, pl);
    call_intersection_global(l, p);
    call_intersection_global(l, r);
    call_intersection_global(l, s);
    call_intersection_global(l, tet);
    call_intersection_global(l, tr);

    call_intersection_global(pl, bbox);
    call_intersection_global(pl, cub);
    call_intersection_global(pl, l);
    call_intersection_global(pl, pl);
    call_intersection_global(pl, p);
    call_intersection_global(pl, r);
    call_intersection_global(pl, s);
    call_intersection_global(pl, sph);
    call_intersection_global(pl, tr);
    call_intersection_global(pl, tet);

    auto plplpl = CGAL::intersection(pl, pl, pl);
    CGAL_USE(plplpl);

    call_intersection_global(p, bbox);
    call_intersection_global(p, cub);
    call_intersection_global(p, l);
    call_intersection_global(p, pl);
    call_intersection_global(p, p);
    call_intersection_global(p, r);
    call_intersection_global(p, s);
    call_intersection_global(p, sph);
    call_intersection_global(p, tr);
    call_intersection_global(p, tet);

    call_intersection_global(r, bbox);
    call_intersection_global(r, cub);
    call_intersection_global(r, l);
    call_intersection_global(r, pl);
    call_intersection_global(r, p);
    call_intersection_global(r, r);
    call_intersection_global(r, s);
    call_intersection_global(r, tr);
    call_intersection_global(r, tet);

    call_intersection_global(s, bbox);
    call_intersection_global(s, cub);
    call_intersection_global(s, l);
    call_intersection_global(s, pl);
    call_intersection_global(s, p);
    call_intersection_global(s, r);
    call_intersection_global(s, s);
    call_intersection_global(s, tet);
    call_intersection_global(s, tr);

    call_intersection_global(sph, pl);
    call_intersection_global(sph, p);
    call_intersection_global(sph, sph);

    call_intersection_global(tr, bbox);
    call_intersection_global(tr, cub);
    call_intersection_global(tr, l);
    call_intersection_global(tr, pl);
    call_intersection_global(tr, p);
    call_intersection_global(tr, r);
    call_intersection_global(tr, s);
    call_intersection_global(tr, tr);
    call_intersection_global(tr, tet);

    call_intersection_global(tet, l);
    call_intersection_global(tet, pl);
    call_intersection_global(tet, p);
    call_intersection_global(tet, r);
    call_intersection_global(tet, s);
    call_intersection_global(tet, tr);

    const auto bbbb = CGAL::intersection(bbox, bbox);
    CGAL_USE(bbbb);

    call_intersection_global(bbox, cub);
    call_intersection_global(bbox, l);
    call_intersection_global(bbox, pl);
    call_intersection_global(bbox, p);
    call_intersection_global(bbox, r);
    call_intersection_global(bbox, s);
    call_intersection_global(bbox, tr);

    call_intersection_global(tet, l);

    // with kernel
    call_intersection_with_kernel(cub, bbox, K());
    call_intersection_with_kernel(cub, cub, K());
    call_intersection_with_kernel(cub, l, K());
    call_intersection_with_kernel(cub, pl, K());
    call_intersection_with_kernel(cub, p, K());
    call_intersection_with_kernel(cub, r, K());
    call_intersection_with_kernel(cub, s, K());
    call_intersection_with_kernel(cub, tr, K());

    call_intersection_with_kernel(l, bbox, K());
    call_intersection_with_kernel(l, cub, K());
    call_intersection_with_kernel(l, l, K());
    call_intersection_with_kernel(l, pl, K());
    call_intersection_with_kernel(l, p, K());
    call_intersection_with_kernel(l, r, K());
    call_intersection_with_kernel(l, s, K());
    call_intersection_with_kernel(l, tr, K());
    call_intersection_with_kernel(l, tet, K());

    call_intersection_with_kernel(pl, bbox, K());
    call_intersection_with_kernel(pl, cub, K());
    call_intersection_with_kernel(pl, l, K());
    call_intersection_with_kernel(pl, pl, K());
    call_intersection_with_kernel(pl, p, K());
    call_intersection_with_kernel(pl, r, K());
    call_intersection_with_kernel(pl, s, K());
    call_intersection_with_kernel(pl, tr, K());
    call_intersection_with_kernel(pl, tet, K());

    //special
    plplpl = K().intersect_3_object()(pl, pl, pl);
    CGAL_USE(plplpl);

    call_intersection_with_kernel(p, bbox, K());
    call_intersection_with_kernel(p, cub, K());
    call_intersection_with_kernel(p, l, K());
    call_intersection_with_kernel(p, pl, K());
    call_intersection_with_kernel(p, p, K());
    call_intersection_with_kernel(p, r, K());
    call_intersection_with_kernel(p, s, K());
    call_intersection_with_kernel(p, sph, K());
    call_intersection_with_kernel(p, tr, K());
    call_intersection_with_kernel(p, tet, K());

    call_intersection_with_kernel(r, bbox, K());
    call_intersection_with_kernel(r, cub, K());
    call_intersection_with_kernel(r, l, K());
    call_intersection_with_kernel(r, pl, K());
    call_intersection_with_kernel(r, p, K());
    call_intersection_with_kernel(r, r, K());
    call_intersection_with_kernel(r, s, K());
    call_intersection_with_kernel(r, tr, K());
    call_intersection_with_kernel(r, tet, K());

    call_intersection_with_kernel(s, bbox, K());
    call_intersection_with_kernel(s, cub, K());
    call_intersection_with_kernel(s, l, K());
    call_intersection_with_kernel(s, pl, K());
    call_intersection_with_kernel(s, p, K());
    call_intersection_with_kernel(s, r, K());
    call_intersection_with_kernel(s, s, K());
    call_intersection_with_kernel(s, tr, K());
    call_intersection_with_kernel(s, tet, K());

    call_intersection_with_kernel(sph, pl, K());
    call_intersection_with_kernel(sph, p, K());
    call_intersection_with_kernel(sph, sph, K());

    call_intersection_with_kernel(tr, bbox, K());
    call_intersection_with_kernel(tr, cub, K());
    call_intersection_with_kernel(tr, l, K());
    call_intersection_with_kernel(tr, pl, K());
    call_intersection_with_kernel(tr, p, K());
    call_intersection_with_kernel(tr, r, K());
    call_intersection_with_kernel(tr, s, K());
    call_intersection_with_kernel(tr, tr, K());
    call_intersection_with_kernel(tr, tet, K());

    call_intersection_with_kernel(tet, l, K());
    call_intersection_with_kernel(tet, pl, K());
    call_intersection_with_kernel(tet, p, K());
    call_intersection_with_kernel(tet, r, K());
    call_intersection_with_kernel(tet, s, K());
    call_intersection_with_kernel(tet, tr, K());

    call_intersection_with_kernel(bbox, cub, K());
    call_intersection_with_kernel(bbox, l, K());
    call_intersection_with_kernel(bbox, pl, K());
    call_intersection_with_kernel(bbox, p, K());
    call_intersection_with_kernel(bbox, r, K());
    call_intersection_with_kernel(bbox, s, K());
    call_intersection_with_kernel(bbox, tr, K());

    // ---------------------------------------------------------------------------------------------
    //                                        DO INTERSECT
    // ---------------------------------------------------------------------------------------------

    call_do_intersect_global(l, cub);
    call_do_intersect_global(l, bbox);
    call_do_intersect_global(l, l);
    call_do_intersect_global(l, pl);
    call_do_intersect_global(l, p);
    call_do_intersect_global(l, r);
    call_do_intersect_global(l, s);
    call_do_intersect_global(l, sph);
    call_do_intersect_global(l, tr);
    call_do_intersect_global(l, tet);

    call_do_intersect_global(pl, bbox);
    call_do_intersect_global(pl, cub);
    call_do_intersect_global(pl, l);
    call_do_intersect_global(pl, pl);
    call_do_intersect_global(pl, p);
    call_do_intersect_global(pl, r);
    call_do_intersect_global(pl, s);
    call_do_intersect_global(pl, sph);
    call_do_intersect_global(pl, tr);
    call_do_intersect_global(pl, tet);

    call_do_intersect_global(p, bbox);
    call_do_intersect_global(p, cub);
    call_do_intersect_global(p, l);
    call_do_intersect_global(p, pl);
    call_do_intersect_global(p, p);
    call_do_intersect_global(p, r);
    call_do_intersect_global(p, s);
    call_do_intersect_global(p, sph);
    call_do_intersect_global(p, tr);
    call_do_intersect_global(p, tet);

    call_do_intersect_global(r, bbox);
    call_do_intersect_global(r, cub);
    call_do_intersect_global(r, l);
    call_do_intersect_global(r, pl);
    call_do_intersect_global(r, p);
    call_do_intersect_global(r, r);
    call_do_intersect_global(r, s);
    call_do_intersect_global(r, tr);
    call_do_intersect_global(r, tet);

    call_do_intersect_global(s, bbox);
    call_do_intersect_global(s, cub);
    call_do_intersect_global(s, pl);
    call_do_intersect_global(s, p);
    call_do_intersect_global(s, l);
    call_do_intersect_global(s, r);
    call_do_intersect_global(s, s);
    call_do_intersect_global(s, sph);
    call_do_intersect_global(s, tr);
    call_do_intersect_global(s, tet);

    call_do_intersect_global(sph, bbox);
    call_do_intersect_global(sph, cub);
    call_do_intersect_global(sph, pl);
    call_do_intersect_global(sph, p);
    call_do_intersect_global(sph, l);
    call_do_intersect_global(sph, r);
    call_do_intersect_global(sph, s);
    call_do_intersect_global(sph, sph);
    call_do_intersect_global(sph, tr);
    call_do_intersect_global(sph, tet);

    call_do_intersect_global(tr, bbox);
    call_do_intersect_global(tr, cub);
    call_do_intersect_global(tr, l);
    call_do_intersect_global(tr, pl);
    call_do_intersect_global(tr, p);
    call_do_intersect_global(tr, r);
    call_do_intersect_global(tr, s);
    call_do_intersect_global(tr, sph);
    call_do_intersect_global(tr, tr);
    call_do_intersect_global(tr, tet);

    call_do_intersect_global(tr, bbox);
    call_do_intersect_global(tr, cub);
    call_do_intersect_global(tr, l);
    call_do_intersect_global(tr, pl);
    call_do_intersect_global(tr, p);
    call_do_intersect_global(tr, r);
    call_do_intersect_global(tr, s);
    call_do_intersect_global(tr, sph);
    call_do_intersect_global(tr, tr);
    call_do_intersect_global(tr, tet);

    call_do_intersect_global(tet, bbox);
    call_do_intersect_global(tet, cub);
    call_do_intersect_global(tet, l);
    call_do_intersect_global(tet, pl);
    call_do_intersect_global(tet, p);
    call_do_intersect_global(tet, r);
    call_do_intersect_global(tet, s);
    call_do_intersect_global(tet, sph);
    call_do_intersect_global(tet, tr);
    call_do_intersect_global(tet, tet);

    call_do_intersect_global(bbox, pl);
    call_do_intersect_global(bbox, l);
    call_do_intersect_global(bbox, r);
    call_do_intersect_global(bbox, s);
    call_do_intersect_global(bbox, tr);
    call_do_intersect_global(bbox, sph);
    call_do_intersect_global(bbox, bbox);

    // with_kernel
    call_do_intersect_with_kernel(l, bbox, K());
    call_do_intersect_with_kernel(l, cub, K());
    call_do_intersect_with_kernel(l, l, K());
    call_do_intersect_with_kernel(l, pl, K());
    call_do_intersect_with_kernel(l, p, K());
    call_do_intersect_with_kernel(l, r, K());
    call_do_intersect_with_kernel(l, s, K());
    call_do_intersect_with_kernel(l, sph, K());
    call_do_intersect_with_kernel(l, tr, K());
    call_do_intersect_with_kernel(l, tet, K());

    call_do_intersect_with_kernel(pl, bbox, K());
    call_do_intersect_with_kernel(pl, cub, K());
    call_do_intersect_with_kernel(pl, l, K());
    call_do_intersect_with_kernel(pl, pl, K());
    call_do_intersect_with_kernel(pl, p, K());
    call_do_intersect_with_kernel(pl, r, K());
    call_do_intersect_with_kernel(pl, s, K());
    call_do_intersect_with_kernel(pl, sph, K());
    call_do_intersect_with_kernel(pl, tr, K());
    call_do_intersect_with_kernel(pl, tet, K());

    call_do_intersect_with_kernel(p, bbox, K());
    call_do_intersect_with_kernel(p, cub, K());
    call_do_intersect_with_kernel(p, l, K());
    call_do_intersect_with_kernel(p, pl, K());
    call_do_intersect_with_kernel(p, p, K());
    call_do_intersect_with_kernel(p, r, K());
    call_do_intersect_with_kernel(p, s, K());
    call_do_intersect_with_kernel(p, sph, K());
    call_do_intersect_with_kernel(p, tr, K());
    call_do_intersect_with_kernel(p, tet, K());

    call_do_intersect_with_kernel(r, bbox, K());
    call_do_intersect_with_kernel(r, cub, K());
    call_do_intersect_with_kernel(r, l, K());
    call_do_intersect_with_kernel(r, pl, K());
    call_do_intersect_with_kernel(r, p, K());
    call_do_intersect_with_kernel(r, r, K());
    call_do_intersect_with_kernel(r, s, K());
    call_do_intersect_with_kernel(r, sph, K());
    call_do_intersect_with_kernel(r, tr, K());
    call_do_intersect_with_kernel(r, tet, K());

    call_do_intersect_with_kernel(s, bbox, K());
    call_do_intersect_with_kernel(s, cub, K());
    call_do_intersect_with_kernel(s, r, K());
    call_do_intersect_with_kernel(s, l, K());
    call_do_intersect_with_kernel(s, p, K());
    call_do_intersect_with_kernel(s, r, K());
    call_do_intersect_with_kernel(s, s, K());
    call_do_intersect_with_kernel(s, sph, K());
    call_do_intersect_with_kernel(s, tr, K());
    call_do_intersect_with_kernel(s, tet, K());

    call_do_intersect_with_kernel(tr, bbox, K());
    call_do_intersect_with_kernel(tr, cub, K());
    call_do_intersect_with_kernel(tr, r, K());
    call_do_intersect_with_kernel(tr, l, K());
    call_do_intersect_with_kernel(tr, p, K());
    call_do_intersect_with_kernel(tr, r, K());
    call_do_intersect_with_kernel(tr, s, K());
    call_do_intersect_with_kernel(tr, tr, K());
    call_do_intersect_with_kernel(tr, tet, K());

    call_do_intersect_with_kernel(tet, bbox, K());
    call_do_intersect_with_kernel(tet, cub, K());
    call_do_intersect_with_kernel(tet, r, K());
    call_do_intersect_with_kernel(tet, l, K());
    call_do_intersect_with_kernel(tet, p, K());
    call_do_intersect_with_kernel(tet, r, K());
    call_do_intersect_with_kernel(tet, s, K());
    call_do_intersect_with_kernel(tet, tr, K());
    call_do_intersect_with_kernel(tet, tet, K());

    call_do_intersect_with_kernel(bbox, cub, K());
    call_do_intersect_with_kernel(bbox, l, K());
    call_do_intersect_with_kernel(bbox, pl, K());
    call_do_intersect_with_kernel(bbox, p, K());
    call_do_intersect_with_kernel(bbox, r, K());
    call_do_intersect_with_kernel(bbox, s, K());
    call_do_intersect_with_kernel(bbox, sph, K());
    call_do_intersect_with_kernel(bbox, tr, K());
  }
}

int main(int argc, char**)
{
  test<CGAL::Cartesian<double> >(argc);
  test<CGAL::Exact_predicates_inexact_constructions_kernel>(argc);
  test<CGAL::Exact_predicates_exact_constructions_kernel>(argc);

  return EXIT_SUCCESS;
}

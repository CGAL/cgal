#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>

#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>

#include <string>
#include <fstream>

typedef const double vec[3];
const vec input[10] =
{ { -0.37503900000000001, 0.5, 0.47710999999999998                     },
  { -0.37503900000000001, 0.5, 0.31415300000000002                     },
  { -0.37503900000000001, 0.5, 0.151196                                },
  { -0.37503900000000001, 0.10979194793565403, 0.151196                },
  { -0.37503899999999996, -0.049577392484619259, 0.151196              },
  { -0.37503900000000001, -0.18584700000000001, 0.151196               },
  { -0.37503900000000001, -0.062128898913781927, 0.20998675245268067   },
  { -0.37503900000000001, 0.11124288530418507, 0.29237289933619037     },
  { -0.37503900000000001, 0.20409758712826165, 0.33649738670770629     },
  { -0.37503900000000001, 0.29695228895233827, 0.38062187407922232     } };

const vec input_normal = { 1., 0., 0. };

template <typename Kernel, typename CDT_2_traits>
bool test(std::string test_name)
{
  std::cerr << "Testing " << test_name << std::endl;
  const unsigned nb_input_points = sizeof(input)/sizeof(vec);

  typedef CGAL::Constrained_triangulation_face_base_2<CDT_2_traits>     CDT_2_fb;
  typedef CGAL::Triangulation_vertex_base_2<CDT_2_traits>               CDT_2_vb;
  typedef CGAL::Triangulation_data_structure_2<CDT_2_vb, CDT_2_fb>      CDT_2_tds;
  typedef CGAL::No_constraint_intersection_requiring_constructions_tag  CDT_2_itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<CDT_2_traits,
                                                     CDT_2_tds,
                                                     CDT_2_itag>        CDT_2;

  typename Kernel::Vector_3 normal(input_normal[0],
                                   input_normal[1],
                                   input_normal[2]);
  CDT_2_traits traits(normal);
  CDT_2 cdt(traits);
  typedef typename CDT_2::Vertex_handle Vh;

  Vh first, previous;
  for(unsigned i = 0; i < nb_input_points; ++i)
  {
    typename Kernel::Point_3 p(input[i][0], input[i][1], input[i][2]);
    Vh current = cdt.insert(p);
    assert(current != previous);
    if(previous != Vh()) cdt.insert_constraint(previous, current);
    if(first == Vh()) first = current;
    previous = current;
  }
  assert(first != Vh());
  assert(previous != first);
  cdt.insert_constraint(previous, first);
  return cdt.is_valid();
}


template <class K, typename CDT_2_traits>
bool test_segment_intersections(std::string test_name)
{
  std::cerr << "Testing " << test_name << std::endl;

  CDT_2_traits traits(typename K::Vector_3(0,0,1));
  typename CDT_2_traits::Intersect_2 intersect = traits.intersect_2_object();

  typedef typename CDT_2_traits::Segment_2 Segment_2;
  typedef typename CDT_2_traits::Point_2 Point_2;

  // test point intersection
  Segment_2 s1( Point_2(0,0,0), Point_2(1,1,0) ),
            s2( Point_2(0,1,0), Point_2(1,0,0) );

  CGAL::Object o = intersect(s1,s2);
  assert( o.is<Point_2>() );

// test segment overlap
  Point_2 pts[4] = { Point_2(0,0,0), Point_2(1,0,1), Point_2(2,0,2), Point_2(4,0,3) };

  o = intersect( Segment_2(pts[0], pts[1]), Segment_2(pts[2], pts[3]) );
  assert( o.empty() );

  // pure overlap
  o = intersect( Segment_2(pts[0], pts[2]), Segment_2(pts[1], pts[3]) );
  assert( o.is<Segment_2>() );
  o = intersect( Segment_2(pts[0], pts[2]), Segment_2(pts[3], pts[1]) );
  assert( o.is<Segment_2>() );
  o = intersect( Segment_2(pts[2], pts[0]), Segment_2(pts[1], pts[3]) );
  assert( o.is<Segment_2>() );
  o = intersect( Segment_2(pts[2], pts[0]), Segment_2(pts[3], pts[1]) );
  assert( o.is<Segment_2>() );
  // segment fully included
  o = intersect( Segment_2(pts[0], pts[3]), Segment_2(pts[1], pts[2]) );
  assert( o.is<Segment_2>() );
  assert( CGAL::object_cast<Segment_2>(o) == Segment_2(pts[1], pts[2]) );
  // segment fully included with shared vertex
  o = intersect( Segment_2(pts[0], pts[1]), Segment_2(pts[0], pts[2]) );
  assert( o.is<Segment_2>() );
  assert( CGAL::object_cast<Segment_2>(o) == Segment_2(pts[0], pts[1]) );
  // segment sharing a vertex
  o = intersect( Segment_2(pts[0], pts[1]), Segment_2(pts[1], pts[2]) );
  assert( o.is<Point_2>() );
  assert( CGAL::object_cast<Point_2>(o) == pts[1]);

  // degenerate segment
  Segment_2 sd(Point_2(1,0,6), Point_2(1,0,7));
  o = intersect( Segment_2(pts[0], pts[2]), sd );
  assert( o.is<Point_2>() );
  o = intersect( Segment_2(pts[0], pts[1]), sd );
  assert( o.is<Point_2>() );
  o = intersect( Segment_2(pts[1], pts[0]), sd );
  assert( o.is<Point_2>() );
  o = intersect( Segment_2(pts[2], pts[3]), sd );
  assert( o.empty() );

  return true;
}

int main()
{
  std::cerr.precision(17);
  typedef CGAL::Exact_predicates_exact_constructions_kernel   Epeck;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;

  bool ok = true;
  ok = ok && test<CGAL::Epick, CGAL::Projection_traits_3<Epick> >("CDT_2 in a 3D plane, with Epick");
  ok = ok && test<CGAL::Epeck, CGAL::Projection_traits_3<Epeck> >("CDT_2 in a 3D plane, with Epeck");
  ok = ok && test_segment_intersections<CGAL::Epick, CGAL::Triangulation_2_projection_traits_3<Epick> >("CDT_2 traits intersection with Epick");
  ok = ok && test_segment_intersections<CGAL::Epeck, CGAL::Triangulation_2_projection_traits_3<Epeck> >("CDT_2 traits intersection with Epeck");

  return ok ? 0 : 1;
}

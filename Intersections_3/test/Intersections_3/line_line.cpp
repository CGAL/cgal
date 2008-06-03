#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Testsuite/assert.h>

typedef CGAL::Quotient< CGAL::MP_Float >                    FT;
typedef CGAL::Cartesian<FT>                                 K;

// The construction test has to be working
// The equal test has to be working

template <class K>
void _test_intersection_construct(K k)
{
  typedef typename K::FT                               FT;
  typedef typename K::Point_3                          Point_3;
  typedef typename K::Line_3                           Line_3;
  typedef typename K::Intersect_3                      Intersect_3;
  typedef typename K::Construct_line_3                 Construct_line_3;

  Intersect_3 theIntersect_3 = k.intersect_3_object();
  Construct_line_3 theConstruct_line_3 = k.construct_line_3_object();

  std::cout << "Testing intersection(Line_arc, Line_arc)..." << std::endl;
  // Testing the case where it overlaps, or do not intersect
  for(int vx=0;vx<4;vx++) {
    for(int vy=1;vy<4;vy++) {
      for(int vz=0;vz<4;vz++) { 
        if(vx == 0 && vy == 0 && vz == 0) continue;
        const FT a = FT(vx);
        const FT b = FT(vy);
        const FT c = FT(vz);
        Line_3 l = theConstruct_line_3(Point_3(0,0,0), Point_3(a,b,c));
        Line_3 l_1 = theConstruct_line_3(Point_3(0,0,0), Point_3(-a,-b,-c));
        Line_3 l_2 = theConstruct_line_3(Point_3(0,0,0), Point_3(-a,-b,7));
        Line_3 l_3 = theConstruct_line_3(Point_3(1,0,0), Point_3(a,b,c));
        Line_3 l_4 = theConstruct_line_3(Point_3(1,0,0), Point_3(a+1,b,c));
        CGAL::Object obj1 = theIntersect_3(l, l_1);
        CGAL::Object obj2 = theIntersect_3(l, l_2);
        CGAL::Object obj3 = theIntersect_3(l, l_3);
        CGAL::Object obj4 = theIntersect_3(l, l_4);
        CGAL::Object obj1l = CGAL::intersection(l, l_1);
        CGAL::Object obj2l = CGAL::intersection(l, l_2);
        CGAL::Object obj3l = CGAL::intersection(l, l_3);
        CGAL::Object obj4l = CGAL::intersection(l, l_4);

        CGAL_test_assert(CGAL::do_intersect(l, l_1));
        CGAL_test_assert(CGAL::do_intersect(l, l_2));
        CGAL_test_assert(CGAL::do_intersect(l, l_3));
        CGAL_test_assert(!CGAL::do_intersect(l, l_4));

        Point_3 interp2, interp3;
        Line_3 interl1;
        CGAL_test_assert(assign(interl1, obj1));
        CGAL_test_assert(assign(interp2, obj2));
        CGAL_test_assert(assign(interp3, obj3));
        CGAL_test_assert(obj4.is_empty());
        CGAL_test_assert(interp2 == Point_3(0,0,0));
        CGAL_test_assert(interp3 == Point_3(a,b,c));
        Point_3 interp5, interp6;
        Line_3 interl2;
        CGAL_test_assert(assign(interl2, obj1l));
        CGAL_test_assert(assign(interp5, obj2l));
        CGAL_test_assert(assign(interp6, obj3l));
        CGAL_test_assert(obj4l.is_empty());
        CGAL_test_assert(interp5 == Point_3(0,0,0));
        CGAL_test_assert(interp6 == Point_3(a,b,c));
      }
    }
  }

  std::cout << "OK!" << std::endl;
}

int main()
{
  K k;
  _test_intersection_construct(k);
  return 0;
}

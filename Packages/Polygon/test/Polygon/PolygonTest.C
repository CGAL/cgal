#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>

typedef CGAL_Cartesian<double> R;
typedef CGAL_Polygon_traits_2<R> Traits;
typedef Traits::Point_2 Point;
typedef Traits::Segment_2 Segment;
typedef Traits::Vector_2 Vector_2;

#include <CGAL/std/list>
#include <CGAL/std/vector>

typedef CGAL_Polygon_2<Traits, CGAL_STD::list<Point> > ListPolygon;
typedef CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > VectorPolygon;

#include <CGAL/std/fstream>
#include <CGAL/std/algorithm>
#include <CGAL/std/cassert>

//-----------------------------------------------------------------------//
//                          test_default_methods
//-----------------------------------------------------------------------//

void test_default_methods(      CGAL_STD::vector<Point>& pvec0,
                          const CGAL_STD::vector<Point>& pvec1,
                                CGAL_STD::list<Point>&   plist0,
                          const CGAL_STD::list<Point>&   plist1)
{
  {
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   x;
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   p0(pvec0.begin(), pvec0.end());
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   p1(pvec0.rbegin(), pvec0.rend());
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   p2(plist0.begin(), plist0.end());
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   p3(plist0.rbegin(), plist0.rend());
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   p0_copy(p0);

    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > y;
    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > p4(pvec0.begin(), pvec0.end());
    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > p5(pvec0.rbegin(), pvec0.rend());
    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > p6(plist0.begin(), plist0.end());
    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > p7(plist0.rbegin(), plist0.rend());
    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > p4_copy(p4);

    x=p0;
    assert(x == p0);
  }

  {
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   x;
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   p0(pvec1.begin(), pvec1.end());
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   p1(pvec1.rbegin(), pvec1.rend());
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   p2(plist1.begin(), plist1.end());
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   p3(plist1.rbegin(), plist1.rend());
    CGAL_Polygon_2<Traits, CGAL_STD::list<Point> >   p0_copy(p0);

    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > y;
    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > p4(pvec1.begin(), pvec1.end());
    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > p5(pvec1.rbegin(), pvec1.rend());
    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > p6(plist1.begin(), plist1.end());
    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > p7(plist1.rbegin(), plist1.rend());
    CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> > p4_copy(p4);

    x=p0;
    assert(x == p0);
  }
}

//-----------------------------------------------------------------------//
//                          test_iterators
//-----------------------------------------------------------------------//

void test_iterators(ListPolygon& p, const ListPolygon& q)
{
  typedef ListPolygon::Vertex_circulator VC;
  typedef ListPolygon::Vertex_const_circulator VCC;
  typedef ListPolygon::Vertex_iterator VI;
  typedef ListPolygon::Vertex_const_iterator VCI;
  typedef ListPolygon::Edge_const_circulator EC;
  typedef ListPolygon::Edge_const_iterator EI;

  CGAL_set_ascii_mode(CGAL_STD::cout);

  {
    VC v = p.vertices_circulator();
    iterator_category(v);
    value_type(v);
    distance_type(v);

    VC vstart(v);
    if (v != 0)
      do {
        CGAL_STD::cout << *v << endl;
        ++v;
      } while (v != vstart);

    for (VI vi = p.vertices_begin(); vi != p.vertices_end(); ++vi)
      CGAL_STD::cout << *vi << endl;

    EC e = p.edges_circulator();
    iterator_category(e);
    value_type(e);
    distance_type(e);

    EC estart(e);
    if (e != 0)
      do {
        CGAL_STD::cout << *e << endl;
        ++e;
      } while (e != estart);

    for (EI ei = p.edges_begin(); ei != p.edges_end(); ++ei)
      CGAL_STD::cout << *ei << endl;
  }

  //-------------------------------------------------------------------//
  {
    VCC v = q.vertices_circulator();
    iterator_category(v);
    value_type(v);
    distance_type(v);

    VCC vstart(v);
    if (v != 0)
      do {
        CGAL_STD::cout << *v << endl;
        ++v;
      } while (v != vstart);

    for (VCI vi = q.vertices_begin(); vi != q.vertices_end(); ++vi)
      CGAL_STD::cout << *vi << endl;

    EC e = q.edges_circulator();
    iterator_category(e);
    value_type(e);
    distance_type(e);

    EC estart(e);
    if (e != 0)
      do {
        CGAL_STD::cout << *e << endl;
        ++e;
      } while (e != estart);

    for (EI ei = q.edges_begin(); ei != q.edges_end(); ++ei)
      CGAL_STD::cout << *ei << endl;
  }
}

//-----------------------------------------------------------------------//
//                          test_stream_operators
//-----------------------------------------------------------------------//

void test_stream_operators(ListPolygon& p)
{
  {
    CGAL_STD::ofstream to("polytest.ascii");
    CGAL_set_ascii_mode(to);
    to << p;
    to.close();

    ListPolygon p_copy;
    CGAL_STD::ifstream from("polytest.ascii");
    CGAL_set_ascii_mode(from);
    from >> p_copy;

    assert(p == p_copy);
  }
  {
    CGAL_STD::ofstream to("polytest.pretty");
    CGAL_set_pretty_mode(to);
    to << p;
  }
  {
    CGAL_STD::ofstream to("polytest.binary");
    CGAL_set_binary_mode(to);
    to << p;
  }
  CGAL_set_pretty_mode(CGAL_STD::cout);
}

//-----------------------------------------------------------------------//
//                          test_access_functions
//-----------------------------------------------------------------------//

void test_access_functions(VectorPolygon& p)
{
  CGAL_STD::cout << "p.size()               = " << p.size() << endl;
  CGAL_STD::cout << "p.is_empty()           = " << (p.is_empty() ? "true" : "false") << endl;

  // test random access methods
#ifndef CGAL_CFG_NO_LAZY_INSTANTIATION
  for (int i=0; i<p.size(); i++) {
    CGAL_STD::cout << "vertex " << i << " = " << p.vertex(i) << endl;
    CGAL_STD::cout << "edge   " << i << " = " << p.edge(i) << endl;
  }

  typedef CGAL_Polygon_2<Traits, CGAL_STD::vector<Point> >::Edge_const_iterator EI;
  EI edges_begin = p.edges_begin();
  EI edges_end   = p.edges_end();
  assert(edges_begin < edges_end);
#endif // CGAL_CFG_NO_LAZY_INSTANTIATION
}

//-----------------------------------------------------------------------//
//                          test_geometric_predicates
//-----------------------------------------------------------------------//

void test_geometric_predicates(ListPolygon& p)
{
  Point point(1.8, 3.133);

  CGAL_STD::cout << "p.is_simple()          = " << int(p.is_simple()) << endl;
  CGAL_STD::cout << "p.is_convex()          = " << int(p.is_convex()) << endl;

  CGAL_Orientation orientation = p.orientation();
  CGAL_STD::cout << "p.orientation()        = ";
  switch (orientation) {
     case CGAL_CLOCKWISE:        CGAL_STD::cout << "clockwise"; break;
     case CGAL_COUNTERCLOCKWISE: CGAL_STD::cout << "counter clockwise"; break;
     case CGAL_COLLINEAR:        CGAL_STD::cout << "collinear"; break;
  }
  CGAL_STD::cout << endl;

  CGAL_Oriented_side oside  = p.oriented_side(point);
  CGAL_STD::cout << "p.oriented_side(point) = ";
  switch (oside) {
    case CGAL_ON_NEGATIVE_SIDE:     CGAL_STD::cout << "on negative side" << endl; break;
    case CGAL_ON_ORIENTED_BOUNDARY: CGAL_STD::cout << "on oriented boundary" << endl; break;
    case CGAL_ON_POSITIVE_SIDE:     CGAL_STD::cout << "on positive side" << endl; break;
  }

  CGAL_STD::cout << "p.bounded_side(point)  = ";
  CGAL_Bounded_side bside   = p.bounded_side(point);
  switch (bside) {
    case CGAL_ON_BOUNDED_SIDE:   CGAL_STD::cout << "on bounded side" << endl; break;
    case CGAL_ON_BOUNDARY:       CGAL_STD::cout << "on boundary" << endl; break;
    case CGAL_ON_UNBOUNDED_SIDE: CGAL_STD::cout << "on unbounded side" << endl; break;
  }

  CGAL_STD::cout << "p.bbox()               = " << p.bbox() << endl;
  CGAL_STD::cout << "p.area()               = " << p.area() << endl;
  CGAL_STD::cout << "*p.left_vertex()       = " << *p.left_vertex() << endl;
  CGAL_STD::cout << "*p.right_vertex()      = " << *p.right_vertex() << endl;
  CGAL_STD::cout << "*p.top_vertex()        = " << *p.top_vertex() << endl;
  CGAL_STD::cout << "*p.bottom_vertex()     = " << *p.bottom_vertex() << endl;
}

//-----------------------------------------------------------------------//
//                          test_update_operations
//-----------------------------------------------------------------------//

void test_update_operations(const ListPolygon& p,
                            const CGAL_STD::vector<Point>& pvec)
{
  // test update functions
  ListPolygon q = p;
  q.reverse_orientation();
  CGAL_STD::cout << "p after reversing orientation: " << q << endl;

  assert(p==p);
  assert(!(p==q));

  typedef ListPolygon::Vertex_iterator VI;
  q=p;
  VI first  = q.vertices_begin();
  VI middle = q.vertices_begin();
  VI last   = q.vertices_end();
  ++middle;
  rotate(first, middle, last);
  assert(p==q);

  // test update operations
  q.push_back(Point(2,3));

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  q.insert(q.vertices_begin(), pvec.begin() + 3, pvec.begin() + 7);
#endif // CGAL_CFG_NO_MEMBER_TEMPLATES
}

//-----------------------------------------------------------------------//
//                          main
//-----------------------------------------------------------------------//

int main()
{
  CGAL_STD::vector<Point> pvec;
  pvec.insert(pvec.end(), Point(0, 3));
  pvec.insert(pvec.end(), Point(1, 3));
  pvec.insert(pvec.end(), Point(1, 4));
  pvec.insert(pvec.end(), Point(3, 4));
  pvec.insert(pvec.end(), Point(3, 2));
  pvec.insert(pvec.end(), Point(2, 3.5));
  pvec.insert(pvec.end(), Point(2, 3));
  pvec.insert(pvec.end(), Point(1.5, 3.5));
  pvec.insert(pvec.end(), Point(1, 2));
  pvec.insert(pvec.end(), Point(4, 1));
  pvec.insert(pvec.end(), Point(7, 3));
  pvec.insert(pvec.end(), Point(6, 5));
  pvec.insert(pvec.end(), Point(4, 2));
  pvec.insert(pvec.end(), Point(4.5, 2));
  pvec.insert(pvec.end(), Point(4, 1.5));
  pvec.insert(pvec.end(), Point(3.5, 3));
  pvec.insert(pvec.end(), Point(3.5, 4.5));
  pvec.insert(pvec.end(), Point(4, 5));
  pvec.insert(pvec.end(), Point(0, 5));

  CGAL_STD::list<Point> plist(pvec.begin(), pvec.end());
  ListPolygon p(pvec.begin(), pvec.end());
  VectorPolygon q(pvec.begin(), pvec.end());

  test_default_methods(pvec, pvec, plist, plist);
  test_iterators(p,p);
  test_stream_operators(p);
  test_access_functions(q);
  test_geometric_predicates(p);
  test_update_operations(p, pvec);

  return 0;
}


#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment;

#include <list>
#include <vector>

using std::vector;
using std::list;
using std::cout;
using std::endl;

typedef CGAL::Polygon_2<K, list<Point> > ListPolygon;
typedef CGAL::Polygon_2<K, vector<Point> > VectorPolygon;

#include <CGAL/Circulator/Circulator_concepts.h>

BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<ListPolygon::Edge_const_circulator>));
BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<ListPolygon::Vertex_circulator>));

#include <fstream>
#include <algorithm>
#include <cassert>

//-----------------------------------------------------------------------//
//                          test_default_methods
//-----------------------------------------------------------------------//

void test_default_methods(      vector<Point>& pvec0,
                          const vector<Point>& pvec1,
                                list<Point>&   plist0,
                          const list<Point>&   plist1)
{
  {
    CGAL::Polygon_2<K, list<Point> >   x;
    CGAL::Polygon_2<K, list<Point> >   p0(pvec0.begin(), pvec0.end());
    // CGAL::Polygon_2<K, list<Point> >   p1(pvec0.rbegin(), pvec0.rend());
    CGAL::Polygon_2<K, list<Point> >   p2(plist0.begin(), plist0.end());
    // CGAL::Polygon_2<K, list<Point> >   p3(plist0.rbegin(), plist0.rend());
    CGAL::Polygon_2<K, list<Point> >   p0_copy(p0);

    CGAL::Polygon_2<K, vector<Point> > y;
    CGAL::Polygon_2<K, vector<Point> > p4(pvec0.begin(), pvec0.end());
    // CGAL::Polygon_2<K, vector<Point> > p5(pvec0.rbegin(), pvec0.rend());
    // CGAL::Polygon_2<K, vector<Point> > p6(plist0.begin(), plist0.end());
    // CGAL::Polygon_2<K, vector<Point> > p7(plist0.rbegin(), plist0.rend());
    CGAL::Polygon_2<K, vector<Point> > p4_copy(p4);

    x=p0;
    assert(x == p0);
  }

  {
    CGAL::Polygon_2<K, list<Point> >   x;
    CGAL::Polygon_2<K, list<Point> >   p0(pvec1.begin(), pvec1.end());
    // CGAL::Polygon_2<K, list<Point> >   p1(pvec1.rbegin(), pvec1.rend());
    CGAL::Polygon_2<K, list<Point> >   p2(plist1.begin(), plist1.end());
    // CGAL::Polygon_2<K, list<Point> >   p3(plist1.rbegin(), plist1.rend());
    CGAL::Polygon_2<K, list<Point> >   p0_copy(p0);

    CGAL::Polygon_2<K, vector<Point> > y;
    CGAL::Polygon_2<K, vector<Point> > p4(pvec1.begin(), pvec1.end());
    // CGAL::Polygon_2<K, vector<Point> > p5(pvec1.rbegin(), pvec1.rend());
    // CGAL::Polygon_2<K, vector<Point> > p6(plist1.begin(), plist1.end());
    // CGAL::Polygon_2<K, vector<Point> > p7(plist1.rbegin(), plist1.rend());
    CGAL::Polygon_2<K, vector<Point> > p4_copy(p4);

    x=p0;
    assert(x == p0);
  }
}

//-----------------------------------------------------------------------//
//                          test_iterators
//-----------------------------------------------------------------------//

void is_input_iterator(std::input_iterator_tag )
{
    return ;
}

void test_iterators(ListPolygon& p, const ListPolygon& q)
{
  typedef ListPolygon::Vertex_circulator VC;
  typedef ListPolygon::Vertex_const_circulator VCC;
  typedef ListPolygon::Vertex_iterator VI;
  typedef ListPolygon::Vertex_const_iterator VCI;
  typedef ListPolygon::Edge_const_circulator EC;
  typedef ListPolygon::Edge_const_iterator EI;

  CGAL::set_ascii_mode(cout);

  {
    VC v = p.vertices_circulator();
    std::iterator_traits<VC>::iterator_category ic1;
    is_input_iterator(ic1);

    VC vstart(v);
    if (v != 0)
      do {
        cout << *v << endl;
        ++v;
      } while (v != vstart);

    for (VI vi = p.vertices_begin(); vi != p.vertices_end(); ++vi)
      cout << *vi << endl;

    EC e = p.edges_circulator();
    std::iterator_traits<VC>::iterator_category ic2;
    is_input_iterator(ic2);

    EC estart(e);
    if (e != 0)
      do {
        cout << *e << endl;
        ++e;
      } while (e != estart);

    for (EI ei = p.edges_begin(); !(p.edges_end() == ei); ++ei) {
      cout << *ei << endl;
      cout << ei->source() << endl;
    }
  }

  //-------------------------------------------------------------------//
  {
    VCC v = q.vertices_circulator();
    std::iterator_traits<VC>::iterator_category ic3;
    is_input_iterator(ic3);

    VCC vstart(v);
    if (v != 0)
      do {
        cout << *v << endl;
        ++v;
      } while (v != vstart);

    for (VCI vi = q.vertices_begin(); vi != q.vertices_end(); ++vi)
      cout << *vi << endl;

    EC e = q.edges_circulator();
    std::iterator_traits<VC>::iterator_category ic4;
    is_input_iterator(ic4);

    EC estart(e);
    if (e != 0)
      do {
        cout << *e << endl;
        ++e;
      } while (e != estart);

    for (EI ei = q.edges_begin(); !(ei == q.edges_end()); ++ei)
      cout << *ei << endl;
  }
}

//-----------------------------------------------------------------------//
//                          test_stream_operators
//-----------------------------------------------------------------------//

void test_stream_operators(ListPolygon& p)
{
  {
    std::ofstream to("polytest.ascii");
    CGAL::set_ascii_mode(to);
    to << p;
    to.close();

    ListPolygon p_copy;
    std::ifstream from("polytest.ascii");
    CGAL::set_ascii_mode(from);
    from >> p_copy;

    assert(p == p_copy);
  }
  {
    std::ofstream to("polytest.pretty");
    CGAL::set_pretty_mode(to);
    to << p;
  }
  {
    std::ofstream to("polytest.binary");
    CGAL::set_binary_mode(to);
    to << p;
  }
  CGAL::set_pretty_mode(cout);
}

//-----------------------------------------------------------------------//
//                          test_access_functions
//-----------------------------------------------------------------------//

void test_access_functions(VectorPolygon& p)
{
  cout << "p.size()               = " << p.size() << endl;
  cout << "p.is_empty()           = " << (p.is_empty() ? "true" : "false") << endl;

  // test random access methods
  for (std::size_t i=0; i<p.size(); i++) {
    cout << "vertex " << i << " = " << p.vertex(i) << endl;
    cout << "edge   " << i << " = " << p.edge(i) << endl;
  }

  typedef CGAL::Polygon_2<K, vector<Point> >::Edge_const_iterator EI;
  EI edges_begin = p.edges_begin();
  EI edges_end   = p.edges_end();
  assert(edges_begin < edges_end);
}

//-----------------------------------------------------------------------//
//                          test_geometric_predicates
//-----------------------------------------------------------------------//

void test_geometric_predicates(ListPolygon& p)
{
  Point point(1.8, 3.133);

  cout << "p.is_simple()          = " << int(p.is_simple()) << endl;
  cout << "p.is_convex()          = " << int(p.is_convex()) << endl;

  CGAL::Orientation orientation = p.orientation();
  cout << "p.orientation()        = ";
  switch (orientation) {
     case CGAL::CLOCKWISE:        cout << "clockwise"; break;
     case CGAL::COUNTERCLOCKWISE: cout << "counter clockwise"; break;
     case CGAL::COLLINEAR:        cout << "collinear"; break;
  }
  cout << endl;

  CGAL::Oriented_side oside  = p.oriented_side(point);
  cout << "p.oriented_side(point) = ";
  switch (oside) {
    case CGAL::ON_NEGATIVE_SIDE:     cout << "on negative side" << endl; break;
    case CGAL::ON_ORIENTED_BOUNDARY: cout << "on oriented boundary" << endl; break;
    case CGAL::ON_POSITIVE_SIDE:     cout << "on positive side" << endl; break;
  }

  cout << "p.bounded_side(point)  = ";
  CGAL::Bounded_side bside   = p.bounded_side(point);
  switch (bside) {
    case CGAL::ON_BOUNDED_SIDE:   cout << "on bounded side" << endl; break;
    case CGAL::ON_BOUNDARY:       cout << "on boundary" << endl; break;
    case CGAL::ON_UNBOUNDED_SIDE: cout << "on unbounded side" << endl; break;
  }

  cout << "p.bbox()               = " << p.bbox() << endl;
  cout << "p.area()               = " << p.area() << endl;
  cout << "*p.left_vertex()       = " << *p.left_vertex() << endl;
  cout << "*p.right_vertex()      = " << *p.right_vertex() << endl;
  cout << "*p.top_vertex()        = " << *p.top_vertex() << endl;
  cout << "*p.bottom_vertex()     = " << *p.bottom_vertex() << endl;
}

//-----------------------------------------------------------------------//
//                          test_update_operations
//-----------------------------------------------------------------------//

void test_update_operations(const ListPolygon& p,
                            const vector<Point>& pvec)
{
  // test update functions
  ListPolygon q = p;
  VectorPolygon pgn(p.vertices_begin(), p.vertices_end());
  q.reverse_orientation();
  cout << "p after reversing orientation: " << q << endl;

  assert(p==p);
  assert(!(p==q));

  typedef ListPolygon::Vertex_iterator VI;
  typedef ListPolygon::Vertex_circulator VC;
  q=p;
  VI middle = q.vertices_begin();
  ++middle;
  q.set(middle, *middle);

  // test update operations
  q.push_back(Point(2,3));
  q.push_back(Point(middle->x(), middle->y()));

  VC c = q.vertices_circulator();
  q.set(c, *middle);
  q.insert(c, Point(2,3)); 
  q.erase(q.vertices_circulator());

  pgn.push_back(Point(pgn.vertices_begin()->x(), 3));
  pgn.set(pgn.vertices_begin(), Point(pgn.vertices_begin()->x(), 3));

  q.insert(q.vertices_begin(), pvec.begin() + 3, pvec.begin() + 7);

  q.insert(q.vertices_circulator(), pvec.begin() + 3, pvec.begin() + 7);

  q.clear();
}

//-----------------------------------------------------------------------//
//                          main
//-----------------------------------------------------------------------//

int main()
{
  vector<Point> pvec;
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

  list<Point> plist(pvec.begin(), pvec.end());
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


//-----------------------------------------------------------------------//
// This example demonstrates how to define a polygon with a custom point
// type with additional attributes.
//
// For this the following is needed:
//
// 1) define a point, vector and segment type
// 2) define a polygon traits class using these types
//
// The bbox() method is only required for compilers that try to instantiate
// all template funcions (g++ 2.7.2 and Borland C++ 5.02).
//-----------------------------------------------------------------------//

#include <CGAL/basic.h>
#include <CGAL/Polygon_2.h>

//#ifdef CGAL_CFG_NO_LAZY_INSTANTIATION
//#include <CGAL/Bbox_2.h>
//#endif // CGAL_CFG_NO_LAZY_INSTANTIATION

//-----------------------------------------------------------------------//
//                          MyPoint
//-----------------------------------------------------------------------//

class MyPoint
{
  private:
    double d_x, d_y;

  public:
    typedef double FT;
    MyPoint(): d_x(0), d_y(0) {}
    MyPoint(double x, double y): d_x(x), d_y(y) {}
    double x() const { return d_x; }
    double y() const { return d_y; }

//#ifdef CGAL_CFG_NO_LAZY_INSTANTIATION
//    CGAL::Bbox_2 bbox() const
//      { return CGAL::Bbox_2(d_x, d_y, d_x, d_y); }
//#endif // CGAL_CFG_NO_LAZY_INSTANTIATION
};

std::ostream& operator<<(std::ostream& to, const MyPoint& p)
{
  return to << "(" << p.x() << ", " << p.y() << ")";
}

//-----------------------------------------------------------------------//
//                          MyVector
//-----------------------------------------------------------------------//

typedef MyPoint MyVector;

MyVector operator-(const MyPoint& p, const MyPoint& q)
{
  return MyVector(p.x() - q.x(), p.y() - q.y());
}

//-----------------------------------------------------------------------//
//                          MySegment
//-----------------------------------------------------------------------//

typedef std::pair<MyPoint,MyPoint> MySegment;

std::ostream& operator<<(std::ostream& to, const MySegment& e)
{
  return to << "source = " << e.first << " target = " << e.second;
}

//-----------------------------------------------------------------------//
//                          MyTraits
//-----------------------------------------------------------------------//

class MyTraits
{
  public:
    typedef double      FT; // the type of the coordinates
    typedef MyPoint     Point_2;
    typedef MyVector    Vector_2;
    typedef MySegment   Segment_2;

    class Less_xy
    {
      public:
        bool operator()(const Point_2& p, const Point_2& q)
        {
          return p.x() < q.x() || (p.x() == q.x() && p.y() < q.y());
        }
    };

    class Less_yx
    {
      public:
        bool operator()(const Point_2& p, const Point_2& q)
        {
          return p.y() < q.y() || (p.y() == q.y() && p.x() < q.x());
        }
    };

    bool lexicographically_xy_smaller(const Point_2& p, const Point_2& q) const
    {
      return p.x() < q.x() || (p.x() == q.x() && p.y() < q.y());
    }

    bool lexicographically_yx_smaller_or_equal(const Point_2& p,
                                               const Point_2& q) const
    {
      return p.y() < q.y() || (p.y() == q.y() && p.x() <= q.x());
    }

    FT cross_product_2(const Vector_2& p, const Vector_2& q) const
    {
      return p.x() * q.y() - q.x() * p.y();
    }

    FT determinant_2(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      return cross_product_2(p-q, p-r);
    }

    int sign(const FT& x) const
    {
      return (x < FT(0)) ? -1: (FT(0) < x) ? 1 : 0;
    }

    CGAL::Orientation orientation(const Point_2& p,
                                 const Point_2& q,
                                 const Point_2& r) const
    {
       int i = sign(determinant_2(p,q,r));
       switch (i) {
       case -1: return CGAL::COUNTERCLOCKWISE;
       case  0: return CGAL::COLLINEAR;
       case  1: return CGAL::CLOCKWISE;
       }

       return CGAL::COLLINEAR; // to prevent compiler warnings
    }

    CGAL::Comparison_result compare_x(const Point_2 &p, const Point_2 &q) const
    {
       if (p.x() == q.x())
          return CGAL::EQUAL;
       return (p.x() < q.x()) ? CGAL::SMALLER : CGAL::LARGER;
    }

    CGAL::Comparison_result compare_y(const Point_2 &p, const Point_2 &q) const
    {
       if (p.y() == q.y())
          return CGAL::EQUAL;
       return (p.y() < q.y()) ? CGAL::SMALLER : CGAL::LARGER;
    }

    bool have_equal_direction(const Vector_2&,
                              const Vector_2&) const
    {
       // N.B. testing for equal directions doesn't make sense when an inexact
       //      numbertype like double is used, so we simply return false
       return false;
    }

    bool do_intersect(const Point_2& p1,
                      const Point_2& q1,
                      const Point_2& p2,
                      const Point_2& q2) const
    {
      // N.B. this implementation is not robust!
      return (sign(determinant_2(p1,q1,p2)) != sign(determinant_2(p1,q1,q2))) &&
             (sign(determinant_2(p2,q2,p1)) != sign(determinant_2(p2,q2,q1)));
    }

    bool is_negative(const FT& x) const
    {
      return (x<0);
    }
};

#include <list>

typedef CGAL::Polygon_2<MyTraits, std::list<MyPoint> > Polygon;
typedef Polygon::Vertex_iterator     VI;
typedef Polygon::Edge_const_iterator EI;

//-----------------------------------------------------------------------//
//                          main
//-----------------------------------------------------------------------//

int main()
{
  // create a polygon and put some points in it
  Polygon p;
  p.push_back(MyPoint(0,0));
  p.push_back(MyPoint(4,0));
  p.push_back(MyPoint(4,4));
  p.push_back(MyPoint(2,2));
  p.push_back(MyPoint(0,4));

  // determine some properties of the polygon
  bool IsSimple    = p.is_simple();
  bool IsConvex    = p.is_convex();
  bool IsClockwise = (p.orientation() == CGAL::CLOCKWISE);
  double Area      = p.area();

  std::cout << "polygon p is";
  if (!IsSimple) std::cout << " not";
  std::cout << " simple." << std::endl;

  std::cout << "polygon p is";
  if (!IsConvex) std::cout << " not";
  std::cout << " convex." << std::endl;

  std::cout << "polygon p is";
  if (!IsClockwise) std::cout << " not";
  std::cout << " clockwise oriented." << std::endl;

  std::cout << "the area of polygon p is " << Area << std::endl;
  std::cout << std::endl;

  // apply some algorithms
  MyPoint q(1,1);
  std::cout << std::endl;

  bool IsInside = (p.bounded_side(q) == CGAL::ON_BOUNDED_SIDE);
  std::cout << "point q is";
  if (!IsInside) std::cout << " not";
  std::cout << " inside polygon p." << std::endl;
  std::cout << std::endl;

  // traverse the vertices and the edges
  int n=0;
  for (VI vi = p.vertices_begin(); vi != p.vertices_end(); ++vi)
    std::cout << "vertex " << n++ << " = " << *vi << std::endl;
  std::cout << std::endl;

  n=0;
  for (EI ei = p.edges_begin(); ei != p.edges_end(); ++ei)
    std::cout << "edge " << n++ << " = " << *ei << std::endl;

  return 0;
}


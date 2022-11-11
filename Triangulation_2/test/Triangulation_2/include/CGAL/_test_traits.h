// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// source        :
// file          : include/CGAL/_test_types.h
// revision      :
// revision_date :
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <CGAL/determinant.h>
#include <CGAL/_test_types.h>
#include <CGAL/enum.h>

namespace CGAL {

// Create a mininal traits class
class Triangulation_test_point {
public:
  typedef Triangulation_test_point Point;
  typedef double TESTFT;
protected:
  double _x, _y;
public:
  Triangulation_test_point() {}
  Triangulation_test_point(double x, double y) : _x(x), _y(y) {}
  Triangulation_test_point(double hx, double hy, double hw) :
    _x(hx/hw), _y(hy/hw)
    {}

  TESTFT test_x() const { return _x; }
  TESTFT test_y() const { return _y; }
  bool   compare(const Point &p) const
    { return test_x()==p.test_x()  &&   test_y()==p.test_y(); }
  bool   uncompare(const Point &p) const { return !compare(p); }
  void   test_set(TESTFT x, TESTFT y) { _x=x; _y=y; }
  bool operator==(const Point &p) const {return this->compare(p);}
};

class Triangulation_test_segment {
  public:
    typedef Triangulation_test_point Point;
  protected:
    Point _p, _q;
  public:
  Triangulation_test_segment() {}
    Triangulation_test_segment(const Point &p, const Point &q)
      : _p(p), _q(q) {}

    void test_set(const Point &p, const Point &q) { _p=p; _q=q; }

 };

class Triangulation_test_line {
  public:
    typedef Triangulation_test_point Point;
  protected:
    Point _p, _q;
  public:
    Triangulation_test_line() {}
    Triangulation_test_line(const Point &p, const Point &q)
      : _p(p), _q(q) {}

  Point first_point() const {return _p;}
  Point second_point() const {return _q;}

//     Triangulation_test_direction direction() {
//             return Triangulation_test_direction(_p,_q);
//     }
//   Triangulation_test_line opposite() {
//             return Triangulation_test_line(_q, _p);
//     }
//     void test_set(const Point &p, const Point &q) { _p=p; _q=q; }
};

class Triangulation_test_direction {
public:
  typedef Triangulation_test_point Point;
  typedef Triangulation_test_line  Line;
protected:
    Point _p, _q;
public:
    Triangulation_test_direction() {}
    Triangulation_test_direction(const Point &p, const Point &q)
      : _p(p), _q(q) {}
    Triangulation_test_direction(const Line &l)
      : _p(l.first_point()), _q(l.second_point()) {}
 //    Triangulation_test_direction perpendicular(const CGAL::Orientation &) const {
//             return *this;
//     }
//     void test_set(const Point &p, const Point &q) { _p=p; _q=q; }
};


class Triangulation_test_ray {
  public:
    typedef Triangulation_test_point Point;
    typedef Triangulation_test_direction Direction;
  protected:
    Point _p;
    Direction _d;
  public:
    Triangulation_test_ray() {}
    Triangulation_test_ray(const Point &p, const Direction &d)
      : _p(p), _d(d) {}
    Triangulation_test_ray(const Point &p, const Point &q)
      : _p(p), _d(p,q) {}
};

class Triangulation_test_triangle {
  public:
    typedef Triangulation_test_point Point;
  protected:
    Point _p, _q, _r;
  public:
    Triangulation_test_triangle() {}
    Triangulation_test_triangle(const Point &p, const Point &q, const Point &r)
      : _p(p), _q(q), _r(r) {}

    const Point &p() const { return _p; }
    const Point &q() const { return _q; }
    const Point &r() const { return _r; }

    void test_set(const Point &p, const Point &q, const Point &r) { _p=p; _q=q; _r=r; }
};


class Triangulation_test_Less_x_2{
public:
  typedef Triangulation_test_point  Point;
  typedef bool                      result_type;

  bool operator()( const Point& p, const Point&  q) const
    {
      return (p.test_x() < q.test_x());
    }
};

class Triangulation_test_Less_y_2{
public:
  typedef Triangulation_test_point  Point;
  typedef bool                      result_type;

  bool operator()( const Point& p, const Point&  q) const
    {
      return (p.test_y() < q.test_y());
    }
};


class Triangulation_test_Compare_x_2{
public:
  typedef Triangulation_test_point  Point;
  typedef CGAL::Comparison_result   result_type;

  CGAL::Comparison_result operator()( const Point& p, const Point&  q) const
    {
      if (p.test_x() < q.test_x()) return CGAL::SMALLER;
      else if (p.test_x() > q.test_x()) return CGAL::LARGER;
      else return CGAL::EQUAL;
    }
};


class Triangulation_test_Compare_y_2{
public:
  typedef Triangulation_test_point  Point;
  typedef CGAL::Comparison_result   result_type;

  CGAL::Comparison_result operator()( const Point& p, const Point&  q) const
    {
      if (p.test_y() < q.test_y()) return CGAL::SMALLER;
      else if (p.test_y() > q.test_y()) return CGAL::LARGER;
      else return CGAL::EQUAL;
    }
};

class Triangulation_test_Orientation_2
{
public:
  typedef Triangulation_test_point  Point;
  typedef CGAL::Orientation   result_type;

  CGAL::Orientation
  operator()( const Point& p, const Point&  q, const Point& r) const
     {
       typedef Point::TESTFT RT;
       RT det = (q.test_x()-p.test_x()) * (r.test_y()-p.test_y())
              - (r.test_x()-p.test_x()) * (q.test_y()-p.test_y());
       if ( det < RT(0) ) return CGAL::CLOCKWISE;
       if ( RT(0) < det ) return CGAL::COUNTERCLOCKWISE;
       return CGAL::COLLINEAR;
     }
};


class Triangulation_test_Side_of_oriented_circle_2
{
public:
  typedef Triangulation_test_point     Point;
  typedef CGAL::Orientation            result_type;

  CGAL::Oriented_side operator() (const Point &p,
                                  const Point &q,
                                  const Point &r,
                                  const Point &t) const
    {
      typedef Point::TESTFT RT;

      RT px( p.test_x());
      RT py( p.test_y());
      RT qx( q.test_x());
      RT qy( q.test_y());
      RT rx( r.test_x());
      RT ry( r.test_y());
      RT tx( t.test_x());
      RT ty( t.test_y());

      RT RT0(0);
      RT RT1(1);

      RT det = CGAL::determinant(px, py, px*px + py*py, RT1,
                                       qx, qy, qx*qx + qy*qy, RT1,
                                       rx, ry, rx*rx + ry*ry, RT1,
                                       tx, ty, tx*tx + ty*ty, RT1);

      return (det<RT0) ? CGAL::ON_NEGATIVE_SIDE
        : ((det==RT0) ? CGAL::ON_ORIENTED_BOUNDARY : CGAL::ON_POSITIVE_SIDE);
    }
};

class Triangulation_test_Construct_circumcenter_2
{
public:
  typedef Triangulation_test_point     Point;
  typedef Point                        result_type;

  Point  operator() (const Point &p,
                     const Point &q,
                     const Point &r) const
    {typedef Point::TESTFT FT;

      FT px( p.test_x());
      FT py( p.test_y());
      FT qx( q.test_x());
      FT qy( q.test_y());
      FT rx( r.test_x());
      FT ry( r.test_y());

      FT px_qx( px - qx);
      FT py_qy( py - qy);
      FT qx_rx( qx - rx);
      FT qy_ry( qy - ry);
      FT rx_px( rx - px);
      FT ry_py( ry - py);

      FT px2_py2( px*px + py*py);
      FT qx2_qy2( qx*qx + qy*qy);
      FT rx2_ry2( rx*rx + ry*ry);

      FT num_x( px2_py2*qy_ry + qx2_qy2*ry_py + rx2_ry2*py_qy);
      FT num_y( px2_py2*qx_rx + qx2_qy2*rx_px + rx2_ry2*px_qx);

      FT den_x( ( px*qy_ry + qx*ry_py + rx*py_qy) * FT( 2));
      FT den_y( ( py*qx_rx + qy*rx_px + ry*px_qx) * FT( 2));

      return( Point( num_x/den_x, num_y/den_y));
    }
};

class Triangulation_test_Construct_midpoint
{
public:
  typedef Triangulation_test_point     Point;
  typedef Point                        result_type;


  Point operator()(const Point &p, const Point &q)
    {
      return Point((p.test_x() + q.test_x())/2. ,
                   (p.test_y() + q.test_y())/2. );
    }
};

class Triangulation_test_Construct_bisector_2
{
public:
  typedef Triangulation_test_point     Point;
  typedef Triangulation_test_line      Line;
  typedef Line                         result_type;

  Line operator()(const Point &p, const Point &q)
    {
      Triangulation_test_Construct_midpoint construct_midpoint;
      Point m  = construct_midpoint(p,q);
      Point mp = Point(m.test_x() + p.test_y() - q.test_y(),
                       m.test_y() - p.test_x() + q.test_x());
      return Line(m,mp);
    }
};

class Triangulation_test_Construct_point_2
{
public:
  typedef Triangulation_test_point     Point;
  typedef const Point&                 result_type;

  const Point& operator()(const Point &p) const
    {
      return p;
    }
};

class Triangulation_test_Construct_segment_2
{
public:
  typedef Triangulation_test_point     Point;
  typedef Triangulation_test_segment   Segment;
  typedef Segment                      result_type;

  Segment operator()(const Point &p, const Point &q)
    {
      return Segment(p,q);
    }
};

class Triangulation_test_Construct_triangle_2
{
 public:
  typedef Triangulation_test_point     Point;
  typedef Triangulation_test_triangle  Triangle;
  typedef Triangle                     result_type;

  Triangle  operator()(const Point &p, const Point &q, const Point& r)
    {
      return Triangle(p,q,r);
    }
};

class Triangulation_test_Construct_direction_2
{
 public:
  typedef Triangulation_test_line       Line;
  typedef Triangulation_test_direction  Direction;
  typedef Direction                     result_type;

  Direction operator()(const Line &l) { return Direction(l);}
};

class Triangulation_test_Construct_ray_2
{
  public:
  typedef Triangulation_test_point     Point;
  typedef Triangulation_test_ray       Ray;
  typedef Triangulation_test_direction Direction;
  typedef Ray                          result_type;

  Ray  operator()(const Point &p, const Direction& d) {return Ray(p,d);}
};


class Triangulation_test_Compare_distance_2
{
public:
  typedef Triangulation_test_point              Point;
  typedef Triangulation_test_Compare_distance_2 Compare_distance_2;
  typedef CGAL::Comparison_result               result_type;

  Triangulation_test_Compare_distance_2() {}

  Point::TESTFT sqr(const Point::TESTFT x) const { return x*x; }
  CGAL::Comparison_result
  operator() ( const Point& p, const Point& q, const Point& r) const
    {
      Point::TESTFT dq = sqr(p.test_x()-q.test_x())
                       + sqr(p.test_y()-q.test_y());
      Point::TESTFT dr = sqr(p.test_x()-r.test_x())
                       + sqr(p.test_y()-r.test_y());
      if( dq < dr) return CGAL::SMALLER ;
      else if (dq == dr) return CGAL::EQUAL;
      return CGAL::LARGER;
    }
};

// class Triangulation_test_Construct_direction_of_line_2
// {
// public:
//   Triangulation_test_Construct_direction_of_line_2(){}
//   typedef Triangulation_test_direction Direction;
//   typedef Triangulation_test_line      Line;
//   Direction operator()(Line l)    { return Direction(l);}
// };

class _Triangulation_test_traits {
public:
  typedef Triangulation_test_point     Point_2;
  typedef Triangulation_test_segment   Segment_2;
  typedef Triangulation_test_line      Line_2;
  typedef Triangulation_test_ray       Ray_2;
  typedef Triangulation_test_triangle  Triangle_2;

  typedef Triangulation_test_Less_x_2       Less_x_2;
  typedef Triangulation_test_Less_y_2       Less_y_2;
  typedef Triangulation_test_Compare_x_2    Compare_x_2;
  typedef Triangulation_test_Compare_y_2    Compare_y_2;
  typedef Triangulation_test_Orientation_2  Orientation_2;
  typedef Triangulation_test_Side_of_oriented_circle_2
                                             Side_of_oriented_circle_2;
  typedef Triangulation_test_Construct_circumcenter_2
                                             Construct_circumcenter_2;
  typedef Triangulation_test_Construct_bisector_2
                                             Construct_bisector_2;
  typedef Triangulation_test_Construct_midpoint
                                             Construct_midpoint;
  typedef Triangulation_test_Compare_distance_2
                                             Compare_distance_2;
  typedef Triangulation_test_Construct_point_2     Construct_point_2;
  typedef Triangulation_test_Construct_segment_2   Construct_segment_2;
  typedef Triangulation_test_Construct_triangle_2  Construct_triangle_2;
  typedef Triangulation_test_Construct_direction_2 Construct_direction_2;
  typedef Triangulation_test_Construct_ray_2       Construct_ray_2;


  Less_x_2
  less_x_2_object() const
    { return Less_x_2();}

  Less_y_2
  less_y_2_object() const
    { return Less_y_2();}

  Compare_x_2
  compare_x_2_object() const
    { return Compare_x_2();}

  Compare_y_2
  compare_y_2_object() const
    { return Compare_y_2();}

  Orientation_2
  orientation_2_object() const
    { return Orientation_2();}

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
    {return Side_of_oriented_circle_2();}

  Construct_circumcenter_2
  construct_circumcenter_2_object() const
    { return Construct_circumcenter_2();}

  Construct_bisector_2
  construct_bisector_2_object() const
    {return Construct_bisector_2();}

  Construct_midpoint
  construct_midpoint_object() const
    {return Construct_midpoint();}

  Compare_distance_2
  compare_distance_2_object() const
    {return Compare_distance_2();}

  Construct_point_2
  construct_point_2_object() const
    { return Construct_point_2();}

  Construct_segment_2
  construct_segment_2_object() const
    {return Construct_segment_2();}

  Construct_triangle_2
  construct_triangle_2_object() const
    {return Construct_triangle_2();}

  Construct_direction_2
  construct_direction_2_object() const
    {return Construct_direction_2();}

  Construct_ray_2
  construct_ray_2_object() const
    {return Construct_ray_2();}
};


std::istream &operator>>(std::istream &is, Triangulation_test_point &p)
{
  typedef Triangulation_test_point::TESTFT FT;
  FT x,y;
  is >> x >> y;
  p.test_set(x,y);
  return is;
}

std::istream &operator>>(std::istream &is, Triangulation_test_triangle &t)
{
  typedef Triangulation_test_triangle::Point Point;
  Point p,q,r;
  is >> p >> q >> r ;
  t.test_set(p,q,r);
  return is;
}


std::ostream &operator<<(std::ostream &os, const Triangulation_test_point &p)
{
  return os << p.test_x() << ' ' << p.test_y();
}


} //namespace CGAL

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
// file          : include/CGAL/_test_types.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <CGAL/determinant.h>
#include <CGAL/_test_types.h>

CGAL_BEGIN_NAMESPACE

// Create a mininal traits class
class Triangulation_test_point {
  public:
    typedef Triangulation_test_point Point;
  protected:
    double _x, _y;
  public:
    typedef double TESTFT;
    Triangulation_test_point() {}
    Triangulation_test_point(const Point &p) : _x(p.test_x()), _y(p.test_y()) {}
    Triangulation_test_point(double x, double y) : _x(x), _y(y) {}
    Triangulation_test_point(double hx, double hy, double hw) : _x(hx/hw), _y(hy/hw) {}

    TESTFT test_x() const { return _x; }
    TESTFT test_y() const { return _y; }
    bool   compare(const Point &p) const { return test_x()==p.test_x() && test_y()==p.test_y(); }
    bool   uncompare(const Point &p) const { return !compare(p); }
    Point  &operator=(const Point &p) { _x=p.test_x(); _y=p.test_y(); return *this; }
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

class Triangulation_test_direction {
  public:
    typedef Triangulation_test_point Point;
  protected:
    Point _p, _q;
  public:
    Triangulation_test_direction() {}
    Triangulation_test_direction(const Point &p, const Point &q)
      : _p(p), _q(q) {}

    Triangulation_test_direction perpendicular(const CGAL::Orientation &o) const {
    	return *this;
    }
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

    Triangulation_test_direction direction() {
    	return Triangulation_test_direction(_p,_q);
    }
  Triangulation_test_line opposite() {
    	return Triangulation_test_line(_q, _p);
    }
    void test_set(const Point &p, const Point &q) { _p=p; _q=q; }
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

class Triangulation_test_distance {
  public:
    typedef Triangulation_test_point Point;
  protected:
    Point _p0, _p1, _p2;
    Point::TESTFT sqr(const Point::TESTFT x) const { return x*x; }
  public:
    Triangulation_test_distance() {}
    Triangulation_test_distance(const Point& p0) : _p0(p0) {}
    Triangulation_test_distance(const Point& p0, const Point& p1) : _p0(p0), _p1(p1) {}
    Triangulation_test_distance(const Point& p0, const Point& p1, const Point& p2)
        : _p0(p0), _p1(p1), _p2(p2) {}
    // misses the constructor with (Pt,Pt,Pt, Gt). To be added later

    void set_point(int i, const Point& p)
    {
        switch(i){
        case 0: _p0 = p; break;
        case 1: _p1 = p; break;
        default: _p2 = p;
        }
    }

    Point get_point(int i) const
    {
      switch(i){
      case 0: return _p0;
      case 1: return _p1;
      }
      return _p2;
    }

    CGAL::Comparison_result compare() const
    {
      Point::TESTFT d01 = sqr(_p0.test_x()-_p1.test_x()) + sqr(_p0.test_y()-_p1.test_y());
      Point::TESTFT d02 = sqr(_p0.test_x()-_p2.test_x()) + sqr(_p0.test_y()-_p2.test_y());
      //the followimg generates a bug with g++2.8.1
      //return CGAL::compare( d01,d02);
      if ( d01 < d02) return CGAL::SMALLER;
      else if ( d01 > d02) return CGAL::LARGER;
      else return CGAL::EQUAL;
    }
};


class _Triangulation_test_traits {
public:
  typedef Triangulation_test_point     Point;
  typedef Triangulation_test_segment   Segment;
  typedef Triangulation_test_line      Line;
  typedef Triangulation_test_ray       Ray;
  typedef Triangulation_test_direction Direction;
  typedef Triangulation_test_triangle  Triangle;

  typedef Triangulation_test_distance  Distance;

  _Triangulation_test_traits() {}
  _Triangulation_test_traits(const _Triangulation_test_traits &) {}
  _Triangulation_test_traits &operator=
        (const _Triangulation_test_traits &) { return *this; }

  static
  CGAL::Comparison_result compare_x(const Point &p, const Point &q)
    { 
      // the following generates a bug with g++ 2.8.1
      // return CGAL::compare(p.test_x(),q.test_x()); 
      if (p.test_x() < q.test_x()) return CGAL::SMALLER;
      else if (p.test_x() > q.test_x()) return CGAL::LARGER;
      else return CGAL::EQUAL;
    }

  static
  CGAL::Comparison_result compare_y(const Point &p, const Point &q)
    {  
      // the following generates a bug with g++ 2.8.1
      // return CGAL::compare(p.test_y(), q.test_y()); 
      if (p.test_y() < q.test_y()) return CGAL::SMALLER;
      else if (p.test_y() > q.test_y()) return CGAL::LARGER;
      else return CGAL::EQUAL;
    }


  static
  bool compare(const Point &p, const Point &q)
    { return p.compare(q); }

  static
  CGAL::Orientation orientation(const Point &p, const Point &q, const Point &r)
    {
      typedef Point::TESTFT RT;
      RT det = (q.test_x()-p.test_x()) * (r.test_y()-p.test_y()) - (r.test_x()-p.test_x()) * (q.test_y()-p.test_y());
      if ( det < RT(0) ) return CGAL::CLOCKWISE;
      if ( RT(0) < det ) return CGAL::COUNTERCLOCKWISE;
      return CGAL::COLLINEAR;
    }

  // side_of_oriented_circle
  static
  CGAL::Oriented_side side_of_oriented_circle(const Point &p, const Point &q,
  					     const Point &r, const Point &t)
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

      RT det = CGAL::det4x4_by_formula(px, py, px*px + py*py, RT1,
                                      qx, qy, qx*qx + qy*qy, RT1,
                                      rx, ry, rx*rx + ry*ry, RT1,
                                      tx, ty, tx*tx + ty*ty, RT1);

      return (det<RT0) ? CGAL::ON_NEGATIVE_SIDE
               : ((det==RT0) ? CGAL::ON_ORIENTED_BOUNDARY : CGAL::ON_POSITIVE_SIDE);
 
    }

  // circumcenter
  static
  Point circumcenter(const Point &p, const Point &q, const Point &r)
    {
      typedef Point::TESTFT FT;

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
    
  // perpendicular bisector (clockwise or counterclockwise)
  static
  Line bisector(const Segment &s)
    {
      return Line();
    }
};

istream &operator>>(istream &is, Triangulation_test_point &p)
{
  Triangulation_test_point::TESTFT x,y;
  is >> x >> y;
  p.test_set(x,y);
  return is;
}

istream &operator>>(istream &is, Triangulation_test_triangle &t)
{
  Triangulation_test_triangle::Point p,q,r;
  is >> p >> q >> r ;
  t.test_set(p,q,r);
  return is;
}




ostream &operator<<(ostream &os, const Triangulation_test_point &p)
{
  return os << p.test_x() << ' ' << p.test_y();
}

//ostream &operator<<(ostream &os, const Triangulation_test_segment &s)
//{
//  return os << s.test_source() << ' ' << s.test_target();
//}

//ostream &operator<<(ostream &os, const Triangulation_test_triangle &t)
//{
//  return os << t.p() <<' '<< t.q() <<' '<< t.r() ;
//}

CGAL_END_NAMESPACE

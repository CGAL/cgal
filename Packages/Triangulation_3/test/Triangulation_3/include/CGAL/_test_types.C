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
// release_date  : 23/10/98
// 
// source        : 
// file          : include/CGAL/_test_types.C
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/determinant.h>
#include <CGAL/_test_types.h>

// Create a mininal traits class
class Triangulation_test_point {
public:
  typedef Triangulation_test_point Point;
protected:
  double _x, _y, _z;
public:
  typedef double TESTFT;
  Triangulation_test_point() {}

  Triangulation_test_point(const Point &p) 
    : _x(p.test_x()), _y(p.test_y()), _z(p.test_z()) {}

  Triangulation_test_point(double x, double y, double z) 
    : _x(x), _y(y), _z(z) {}

  Triangulation_test_point(double hx, double hy, double hz, double hw) 
    : _x(hx/hw), _y(hy/hw) , _z(hz/hw) {}

  TESTFT test_x() const { return _x; }
  TESTFT test_y() const { return _y; }
  TESTFT test_z() const { return _z; }
  bool   compare(const Point &p) const 
    { return test_x()==p.test_x() 
	&& test_y()==p.test_y() 
	&& test_z()==p.test_z() ; }

  bool   uncompare(const Point &p) const { return ( ! compare(p)); }
  Point  &operator=(const Point &p) 
    { _x=p.test_x(); _y=p.test_y();_z=p.test_z();  return *this; }
  void   test_set(TESTFT x, TESTFT y, TESTFT z) { _x=x; _y=y; _z=z;}
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
  bool operator==(const Triangulation_test_segment &t) const
    {
      return((t._p==this->_p) && (t._q==this->_q));
    }

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

  Triangulation_test_direction perpendicular(const CGAL::Orientation &) const {
    return *this;
  }
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
  void test_set(Point p, Point q, Point r) {_p=p ; _q=q ; _r=r;}
  bool operator==(const Triangulation_test_triangle &t) const
    {
      return((t._p==this->_p) && (t._q==this->_q)&&(t._r==this->_r));
    }
  
};

class Triangulation_test_tetrahedron {
public:
  typedef Triangulation_test_point Point;
protected:
  Point _p, _q, _r, _s;
public:
  Triangulation_test_tetrahedron() {}
  Triangulation_test_tetrahedron(const Point &p, const Point &q, const Point &r, const Point &s)
    : _p(p), _q(q), _r(r), _s(s) {}

  const Point &p() const { return _p; }
  const Point &q() const { return _q; }
  const Point &r() const { return _r; }
  const Point &s() const { return _s; }

  bool operator==(const Triangulation_test_tetrahedron &t) const
    {
      return((t._p==this->_p) && (t._q==this->_q)&&(t._r==this->_r)&&(t._s==this->_s));
    }
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
  Point::TESTFT d01 = sqr(_p0.test_x()-_p1.test_x()) + sqr(_p0.test_y()-_p1.test_y()) + sqr(_p0.test_z()-_p1.test_z()) ;
  Point::TESTFT d02 = sqr(_p0.test_x()-_p2.test_x()) + sqr(_p0.test_y()-_p2.test_y()) + sqr(_p0.test_z()-_p1.test_z());
  return CGAL_NTS compare( d01,d02);
}
};

class _Triangulation_test_traits_3 {
public :
  typedef Triangulation_test_point                Point;

  typedef CGAL::Point_2<CGAL::Cartesian<double> >        Point2;
  typedef CGAL::Point_3<CGAL::Cartesian<double> >     Point3;

  typedef Triangulation_test_segment              Segment;
  typedef Triangulation_test_line                 Line;
  typedef Triangulation_test_ray                  Ray;
  typedef Triangulation_test_direction            Direction;
  typedef Triangulation_test_triangle             Triangle;
  typedef Triangulation_test_tetrahedron          Tetrahedron;
  typedef Triangulation_test_distance             Distance;
 


  _Triangulation_test_traits_3() {}
  _Triangulation_test_traits_3(const _Triangulation_test_traits_3 &) {}
_Triangulation_test_traits_3 &operator=
(const _Triangulation_test_traits_3 &) { return *this; }

static
CGAL::Comparison_result compare_x(const Point &p, const Point &q)
{ return CGAL_NTS compare(p.test_x(),q.test_x()); }

static
CGAL::Comparison_result compare_y(const Point &p, const Point &q)
{ return CGAL_NTS compare(p.test_y(), q.test_y()); }

static
CGAL::Comparison_result compare_z(const Point &p, const Point &q)
{ return CGAL_NTS compare(p.test_z(), q.test_z()); }

static
bool equal(const Point &p, const Point &q)
{ return p.compare(q); }

static
CGAL::Orientation orientation(const Point & p, const Point & q, const
			      Point & r,const Point & s) 
{
  Point3 p1(p.test_x(),p.test_y(),p.test_z());
  Point3 p2(q.test_x(),q.test_y(),q.test_z());
  Point3 p3(r.test_x(),r.test_y(),r.test_z());
  Point3 p4(s.test_x(),s.test_y(),s.test_z());
  return CGAL::orientation(p1,p2,p3,p4);
}

static      
CGAL::Orientation orientation_in_plane(const Point &q, const Point &r, 
				       const Point &s, const Point &p)
{

  Point2 pxy(p.test_x(), p.test_y());
  Point2 qxy(q.test_x(), q.test_y());
  Point2 rxy(r.test_x(), r.test_y());
  Point2 sxy(s.test_x(), s.test_y());
  CGAL::Orientation oxy_qrs = CGAL::orientation(qxy,rxy,sxy);

  if ( oxy_qrs !=  CGAL::COLLINEAR ) {
    // the projection on x,y is OK
    // tests whether pxy is on the same side of qxy, rxy as sxy
    CGAL::Orientation oxy_qrp = CGAL::orientation(qxy,rxy,pxy);
    if ( oxy_qrp == oxy_qrs) {
      return CGAL::POSITIVE;
    }
    else {
      if ( oxy_qrp == CGAL::COLLINEAR ) { return CGAL::COLLINEAR; }
      else { return CGAL::NEGATIVE; }
    }
  }

  // else : must project on another plane
  // tests on which plane :

  if ( ( qxy.x() != rxy.x() ) || 
       ( qxy.x() != sxy.x() ) ) {
    // projection on x,z-plane is ok
    Point2 pxz(p.test_x(), p.test_z());
    Point2 qxz(q.test_x(), q.test_z());
    Point2 rxz(r.test_x(), r.test_z());
    Point2 sxz(s.test_x(), s.test_z());
    // tests whether pxz is on the same side of qxz, rxz as sxz
    CGAL::Orientation oxz_qrs = CGAL::orientation(qxz,rxz,sxz);
    CGAL::Orientation oxz_qrp = CGAL::orientation(qxz,rxz,pxz);
    if ( oxz_qrp == oxz_qrs) {
      return CGAL::POSITIVE;
    }
    else {
      if ( oxz_qrp == CGAL::COLLINEAR ) { return CGAL::COLLINEAR;	}
      else { return CGAL::NEGATIVE; }
    }
  }
   
  // else : projection on y,z-plane
  Point2 pyz(p.test_y(), p.test_z());
  Point2 qyz(q.test_y(), q.test_z());
  Point2 ryz(r.test_y(), r.test_z());
  Point2 syz(s.test_y(), s.test_z());
  // tests whether pyz is on the same side of qyz, ryz as syz
  CGAL::Orientation oyz_qrs = CGAL::orientation(qyz,ryz,syz);
  CGAL::Orientation oyz_qrp = CGAL::orientation(qyz,ryz,pyz);
  if ( oyz_qrp == oyz_qrs) {
    return CGAL::POSITIVE;
  }
  else {
    if ( oyz_qrp == CGAL::COLLINEAR ) { return CGAL::COLLINEAR;	}
    else { return CGAL::NEGATIVE; }
  }
}

bool collinear(const Point & p,
	       const Point & q,
	       const Point & r) const
{
  Point3 p1(p.test_x(),p.test_y(),p.test_z());
  Point3 p2(q.test_x(),q.test_y(),q.test_z());
  Point3 p3(r.test_x(),r.test_y(),r.test_z());

  return CGAL::collinear(p1,p2,p3);
}
};


std::istream &operator>>(std::istream &is, Triangulation_test_point &p)
{
  Triangulation_test_point::TESTFT x,y, z;
  is >> x >> y >> z;
  p.test_set(x,y,z);
  return is;
}

std::istream &operator>>(std::istream &is, Triangulation_test_triangle &t)
{
  Triangulation_test_triangle::Point p,q,r;
  is >> p >> q >> r ;
  t.test_set(p,q,r);
  return is;
}




std::ostream &operator<<(std::ostream &os, const Triangulation_test_point &p)
{
  return os << p.test_x() << ' ' << p.test_y() << ' ' << p.test_z() ;
}

// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Parabola_2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_PARABOLA_2_H
#define CGAL_PARABOLA_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#include <CGAL/Cartesian.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/determinant.h>
#include <CGAL/IO/Window_stream.h>


CGAL_BEGIN_NAMESPACE


template < class Point, class Weight, class Line >
class Parabola_2
{
public:
  typedef CGAL::Weighted_point< Point, Weight >   Weighted_point;
  typedef CGAL::Point_2< Cartesian<double> >      Point_2;
  typedef CGAL::Segment_2< Cartesian<double> >    Segment_2;
  typedef CGAL::Line_2< Cartesian<double> >       Line_2;

protected:
  // static stuff
  static const double STEP;

  inline static
  double square(const double &x)
  {
    return x * x;
  }

  inline static
  double norm2(const Point_2& p)
  {
    return square(p.x()) + square(p.y());
  }

  inline static
  double distance2(const Point_2& p1, const Point_2& p2)
  {
    double dx = p1.x()-p2.x();
    double dy = p1.y()-p2.y();
    return square(dx) + square(dy);
  }

  inline static
  double distance(const Point_2& p1, const Point_2& p2)
  {
    return sqrt( distance2(p1, p2) );
  }

  inline static
  double distance(const Point_2& p, const Line_2& l)
  {
    return ( p.x() * l.a() + p.y() * l.b() + l.c() ) /
      sqrt( square(l.a()) + square(l.b()) );
  }

  // instance stuff
  Point_2 c;
  Line_2 l;
  Point_2 o;

  inline
  Point_2 lchain(const double &t) const
  {
    std::vector< Point_2 > p = compute_points(t);
    if ( right(p[0]) )  return p[1];
    return p[0];
  }

  inline
  Point_2 rchain(const double &t) const
  {
    std::vector< Point_2 > p = compute_points(t);
    if ( right(p[0]) )  return p[0];
    return p[1]; 
  }

  std::vector< Point_2 > compute_points(const double &d) const
  {
    assert(d >= 0);
    double d1 = distance(o, c) + d;
    double d2 = distance(o, l) + d;
    d2 = d1;
    d1 *= d1;

    std::vector< Point_2 > p;

    if ( l.a() == ZERO ) {
      double y = d2 * CGAL_NTS sign(l.b()) - l.c() / l.b();

      double C = CGAL_NTS square(y) - 2 * c.y() * y + 
	square(c.x()) + square(c.y()) - d1;

      double D = square(c.x()) - C;

      D = CGAL_NTS abs(D);

      double x1 = CGAL_NTS sqrt(D) + c.x();
      double x2 = -CGAL_NTS sqrt(D) + c.x();

      p.push_back(Point_2(x1, y));
      p.push_back(Point_2(x2, y));

      return p;
    }

    double A = d2 * sqrt( square(l.a()) + square(l.b()) ) - l.c();
    double B = square(c.x()) + square(c.y()) - d1;

    double alpha = 1 + square(l.b() / l.a());
    double beta = A * l.b() / square(l.a()) + c.y()
      - c.x() * l.b() / l.a();
    double gamma = square(A / l.a()) + B - 2 * c.x() * A / l.a();

    double D = square(beta) - alpha * gamma;

    D = CGAL_NTS abs(D);

    double y1 = (beta + sqrt(D)) / alpha;
    double y2 = (beta - sqrt(D)) / alpha;

    double x1 = (A - l.b() * y1) / l.a();
    double x2 = (A - l.b() * y2) / l.a();

    p.push_back(Point_2(x1, y1));
    p.push_back(Point_2(x2, y2));

    return p;
  }

  bool right(const Point_2& p) const
  {
    return CGAL_NTS
      is_positive( det3x3_by_formula< double >(c.x(), c.y(), 1,
					       o.x(), o.y(), 1,
					       p.x(), p.y(), 1) );
  }

  inline
  Point_2 midpoint(const Point_2& p1, const Point_2& p2) const
  {
    double t1 = t(p1);
    double t2 = t(p2);
    double midt = (t1+t2)/2;
    return f(midt);
  }

  inline
  Point_2 f(double t) const
  {
    if ( CGAL_NTS is_negative(t) )  return rchain(-t);
    return lchain(t);
  }

  inline
  double t(const Point_2 &p) const
  {
    double tt = distance(p, c) - distance(c, o);
    if ( right(p) )  return -tt;
    return tt;
  }

  void compute_origin()
  {
    double d = (l.a() * c.x() + l.b() * c.y() + l.c())
      / (  2 * ( square(l.a()) + square(l.b()) )  );
    o = Point_2(c.x() - l.a() * d, c.y() - l.b() * d);
  }

public:
  Parabola_2() {}

  Parabola_2(const Weighted_point &p, const Line &l)
  {
    this->c = Point_2(CGAL_NTS to_double(p.x()),
		      CGAL_NTS to_double(p.y()));
    double d_a = CGAL_NTS to_double(l.a());
    double d_b = CGAL_NTS to_double(l.b());
    double len = CGAL_NTS sqrt(CGAL_NTS square(d_a) +
			       CGAL_NTS square(d_b));

    double r = CGAL_NTS to_double(p.weight()) * len;

    this->l = Line_2(CGAL_NTS to_double(-l.a()),
		     CGAL_NTS to_double(-l.b()),
		     CGAL_NTS to_double(-l.c()) + r);
    compute_origin();
  }

  Oriented_side
  side_of_parabola(const Point& p) const
  {
    Point_2 q(CGAL_NTS to_double(p.x()),
	      CGAL_NTS to_double(p.y()));

    double d = distance(q, c) - fabs(distance(q, l));
    if ( d < 0 )  return ON_NEGATIVE_SIDE;
    if ( d > 0 )  return ON_POSITIVE_SIDE;
    return ON_ORIENTED_BOUNDARY;
  }


  inline Line_2 line() const
  {
    return l;
  }

  inline Point_2 center() const
  {
    return c;
  }

  template< class Stream >
  void draw(Stream& W) const
  {
    std::vector< Point_2 > p;
    std::vector< Point_2 > pleft, pright;

    pleft.push_back(o);
    pright.push_back(o);

    for (int i = 1; i <= 100; i++) {
      p = compute_points(i * i * STEP);

      W << p[0];
      W << p[1];

      if ( p.size() > 0 ) {
	if ( right(p[0]) ) {
	  pright.push_back(p[0]);
	  pleft.push_back(p[1]);
	} else {
	  pright.push_back(p[1]);
	  pleft.push_back(p[0]);
	}
      }
    }

    for (unsigned int i = 0; i < pleft.size() - 1; i++) {
      W << Segment_2(pleft[i], pleft[i+1]);
    }

    for (unsigned int i = 0; i < pright.size() - 1; i++) {
      W << Segment_2(pright[i], pright[i+1]);
    }

    W << o;
  }
};

template < class Point, class Weight, class Line >
const double Parabola_2<Point,Weight,Line>::STEP = 2;

template< class Stream, class Point, class Weight, class Line >
inline
Stream& operator<<(Stream& s,
		   const Parabola_2< Point, Weight, Line > &P)
{
  P.draw(s);
  return s;
}




CGAL_END_NAMESPACE

#endif // CGAL_PARABOLA_2_H

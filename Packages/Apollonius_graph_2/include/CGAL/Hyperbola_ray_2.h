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
// file          : include/CGAL/Hyperbola_ray_2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_HYPERBOLA_RAY_2_H
#define CGAL_HYPERBOLA_RAY_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/Hyperbola_segment_2.h>

CGAL_BEGIN_NAMESPACE

template < class Point, class Weight >
class Hyperbola_ray_2 : public Hyperbola_segment_2< Point, Weight >
{
public:
  typedef Sign                                   Hyperbola_direction;
  typedef CGAL::Weighted_point< Point, Weight >  Weighted_point;
  //  typedef typename R::RT         FT;
  typedef double                                 FT;
  typedef CGAL::Point_2< Cartesian<double> >     Point_2;
  typedef CGAL::Segment_2< Cartesian<double> >   Segment_2;
  typedef CGAL::Ray_2< Cartesian<double> >       Ray_2;

protected:
  static const FT OFFSET;

  template< class Stream >
  inline
  void draw_ray(Stream &W) const
  {
    W << Ray_2(this->p1, this->p2);
  }

  Weighted_point _f1, _f2;
  Point _p;
  Hyperbola_direction _dir;

public:
  Hyperbola_ray_2() : Hyperbola_segment_2< Point, Weight >() {}


  Hyperbola_ray_2(const Weighted_point &f1,
		  const Weighted_point &f2,
		  const Point &p,
		  const Hyperbola_direction& direction) :
    Hyperbola_segment_2< Point, Weight >(f1, f2, p, p),
    _f1(f1), _f2(f2), _p(p), _dir(direction)
  {
    FT t1 = t(this->p1);
    if ( direction == POSITIVE ) {
      this->p2 = f(t1 + this->STEP * OFFSET);
    } else {
      this->p2 = f(t1 - this->STEP * OFFSET);
    }
  }


#if defined CGAL_QT_WIDGET_H
  inline
  void draw_qt(Qt_widget& s)
  {
    if ( CGAL_NTS is_zero(r) ) {
      draw_ray(s);
      return;
    }

    double width = s.x_max() - s.x_min();
    double height = s.y_max() - s.y_min();

    if ( width > height ) {
      STEP = height / 100.0;
    } else {
      STEP = width / 100.0;
    }

    FT t1 = t(this->p1);
    if ( _dir == POSITIVE ) {
      this->p2 = f(t1 + STEP * OFFSET);
    } else {
      this->p2 = f(t1 - STEP * OFFSET);
    }
    
    Hyperbola_segment_2< Point, Weight >::draw(s);
  }
#endif

  template< class Stream >
  inline
  void draw(Stream& s) const
  {
    if ( CGAL_NTS is_zero(this->r) ) {
      draw_ray(s);
      return;
    }

    Hyperbola_segment_2< Point, Weight >::draw(s);
  }
  
};

template < class Point, class Weight >
const double Hyperbola_ray_2<Point,Weight>::OFFSET = 1000;



template< class Stream, class Point, class Weight >
inline
Stream& operator<<(Stream &s,
		   const Hyperbola_ray_2< Point, Weight > &H)
{
  H.draw(s);
  return s;
}

#if defined CGAL_QT_WIDGET_H
template< class Point, class Weight >
inline
Qt_widget& operator<<(Qt_widget &s,
		      Hyperbola_ray_2< Point, Weight > &H)
{
  H.draw_qt(s);
  return s;
}
#endif

CGAL_END_NAMESPACE

#endif // CGAL_HYPERBOLA_RAY_2_H

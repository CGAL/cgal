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
// release_date  : 1999, October 01
//
// file          : include/CGAL/min_sqr_distance_traits.h
// package       : bops (2.2)
// source        : include/CGAL/min_sqr_distance_traits.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_MIN_SQR_DISTANCE_TRAITS_H
#define CGAL_MIN_SQR_DISTANCE_TRAITS_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif
#ifndef CGAL_HOMOGENEOUS_H
#include <CGAL/Homogeneous.h>
#endif
#include <list>
#include <vector>
#include <set>
#include <algorithm>
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif

CGAL_BEGIN_NAMESPACE


  struct _dPoint {
    _dPoint() : _x(0.0), _y(0.0) {}
    _dPoint(double sx, double sy) : _x(sx), _y(sy) {}
    _dPoint(const _dPoint& p) : _x(p.x()), _y(p.y()) {}
    double x() const { return _x; }
    double y() const { return _y; }
    _dPoint& operator=(const _dPoint& p) {
      _x= p.x();
      _y= p.y();
      return *this;
    }
    double _x, _y;
  };


template < class _R >
struct min_sqr_distance_traits : public _R {
  typedef _R R;
  typedef typename _R::FT NT;
  //typedef Point_2< Cartesian<double> > dPoint;
  typedef _dPoint dPoint;

  typedef CGAL::Point_2<R>  Point;


  NT sqr( const NT& a ) const { return a * a; }
  double dsqr( const double& a ) const { return a * a; }

  double to_double(const NT& a) const {
    return CGAL::to_double(a);
  }
 

  dPoint to_dPoint(const Point& pt) const {
    return dPoint(to_double(pt.x()), to_double(pt.y()) );
  }


  double double_squared_distance(const Point& p1, const Point& p2) const {
    /* returns the square go the euclidean metric: 
               (p1.x-p2.x)^2 + (p1.y-p2.y)^2
    */
    return to_double( sqr(p1.x()-p2.x()) + sqr(p1.y()-p2.y()) );
  }

  double squared_distance(const dPoint& p1, const dPoint& p2) const {
    /* returns the square go the euclidean metric: 
               (p1.x-p2.x)^2 + (p1.y-p2.y)^2
    */
    return dsqr(p1.x()-p2.x()) + dsqr(p1.y()-p2.y());
  }

  double squared_distance_x(const dPoint& p1, const dPoint& p2) const {
    /* returns (p1.x-p2.x)^2 */
    return dsqr(p1.x()-p2.x());
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_MIN_SQR_DISTANCE_TRAITS_H

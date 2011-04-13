// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Weighted_point.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//                 Sylvain Pion
//
// coordinator   : Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_WEIGHTED_POINT_H
#define CGAL_WEIGHTED_POINT_H

CGAL_BEGIN_NAMESPACE

template < class Pt, class We >
class Weighted_point : public Pt
{
public:
  typedef We Weight;
  typedef Pt Point;
  typedef typename Point::RT RT;

  Weighted_point (const Point &p=Point(), const Weight &w = Weight(0))
      : Point(p), _weight(w) {}

  Point point() const
  {
      return (Point)*this;
  }

  Weight weight() const
  {
      return _weight;
  }

// The following power() member functions are not used at the moment.
// They belong to the traits class anyway.
//
//  Weight power(const Point &p)
//  {	
//      return squared_distance(*this, p) - weight();
//  }
// 
//  Weight power(const Weighted_point &p)
//  {	
//      return squared_distance(*this, p) - weight() - p.weight();
//  }

private:
  Weight _weight;
};


template < class Point, class Weight >
std::ostream &
operator<<(std::ostream &os, const Weighted_point<Point,Weight> &p)
{
	return os << p.point() << " " << p.weight();
}

template < class Point, class Weight >
std::istream &
operator>>(std::istream &is, Weighted_point<Point,Weight> &wp)
{
	Weight w;
	Point p;
	is >> p >> w;
	wp = Weighted_point<Point,Weight>(p,w);
	return is;
}

CGAL_END_NAMESPACE

#endif // CGAL_WEIGHTED_POINT_H

// Copyright (c) 1999-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Mariette Yvinec
//                 Sylvain Pion

#ifndef CGAL_WEIGHTED_POINT_H
#define CGAL_WEIGHTED_POINT_H

CGAL_BEGIN_NAMESPACE

template < class Pt, class We >
class Weighted_point : public Pt
{
public:
  typedef We Weight;
  typedef Pt Point;

  Weighted_point ()
      : Point(), _weight(0) {}

  //explicit
  Weighted_point (const Point &p)
      : Point(p), _weight(0) {}

  Weighted_point (const Point &p, const Weight &w)
      : Point(p), _weight(w) {}

  const Point & point() const
  {
      return *this;
  }

  const Weight & weight() const
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
	if (is)
	    wp = Weighted_point<Point,Weight>(p,w);
	return is;
}

CGAL_END_NAMESPACE

#endif // CGAL_WEIGHTED_POINT_H

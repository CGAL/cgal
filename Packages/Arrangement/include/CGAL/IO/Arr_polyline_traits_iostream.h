// Copyright (c) 2001  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifdef CGAL_ARR_POLYLINE_TRAITS_H
#ifndef CGAL_ARR_POLYLINE_TRAITS_IOSTREAM_H   
#define CGAL_ARR_POLYLINE_TRAITS_IOSTREAM_H  


CGAL_BEGIN_NAMESPACE

/*!
 * Output operator for a polyline.
 */
template <class Segment_traits_>
::std::ostream& operator<< (::std::ostream& os,
                            const Polyline_2<Segment_traits_>& pl)
{
  typedef Polyline_2<Segment_traits_>  Curve_2;
  
  typename Curve_2::const_iterator it;

  // Print out the number of points in the polyline.
  os << pl.points();

  // Print out the polyline points.
  for (it = pl.begin(); it != pl.end(); it++) 
    os << " " << (*it);

  return (os);
}

/*!
 * Input operator for a polyline.
 */
template <class Segment_traits_>
::std::istream& operator>> (::std::istream& is,
                            Polyline_2<Segment_traits_>& pl)
{
  typedef Polyline_2<Segment_traits_>  Curve_2;
  typedef typename Curve_2::Point_2    Point_2;

  // Read the number of input points.
  int    n_pts;

  is >> n_pts;

  // Read n_pts points to a list.
  Point_2              p;
  ::std::list<Point_2> pts;
  int                  i;

  for (i = 0; i < n_pts; i++)
  {
    is >> p;
    pts.push_back(p);
  }

  // Create the polyline curve.
  pl = Curve_2 (pts.begin(), pts.end());

  return (is);
}

CGAL_END_NAMESPACE

#endif
#endif

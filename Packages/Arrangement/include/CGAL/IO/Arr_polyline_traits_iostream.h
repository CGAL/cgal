// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-62 $
// release_date  : $CGAL_Date: 2001/05/11 $
//
// file          : include/CGAL/IO/Arr_polyline_traits_iostream.h
// package       : Arrangement (1.82)
// maintainer    : Efi Fogel <efif@post.tau.ac.il>
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <danha@post.tau.ac.il>)
//
// ======================================================================

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

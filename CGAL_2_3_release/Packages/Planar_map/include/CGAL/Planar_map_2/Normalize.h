// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : include/CGAL/Planar_map_2/Normalize.h
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ============================================================================

#ifndef CGAL_PLANAR_MAP_2_NORMALIZE_H
#define CGAL_PLANAR_MAP_2_NORMALIZE_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_POINT_2_H
template <class R_> // default does nothing.
inline bool normalize_coordinates(Point_2<R_>& pt){
    return false;
  }
#ifdef CGAL_CARTESIAN_H
#ifdef CGAL_LEDA_RATIONAL_H
// remove inline in future use
  inline bool normalize_coordinates(Point_2<Cartesian<leda_rational> >& pt)  {
    typedef CGAL::Point_2<Cartesian<leda_rational> > Point;
    typedef Point::RT RT;
    RT x, y;
    x = pt.x();
    y = pt.y();
    if (CGAL_NTS is_zero(x) &&  CGAL_NTS is_zero(y)) return false;
    x.normalize();
    y.normalize();
    pt=Point(x,y);
    return true;
  }  
#endif
#ifdef CGAL_LEDA_INTEGER_H
// remove inline in future use
  inline bool normalize_coordinates(Point_2<Cartesian<leda_integer> >& pt)  {
    typedef typename CGAL::Point_2<Cartesian<leda_integer> >::RT RT;
    RT g, x, y;
    x = pt.x();
    y = pt.y();
    if (CGAL_NTS is_zero(x) &&  CGAL_NTS is_zero(y)) return false;
    else {
      g = gcd(x,y);
      if (g != RT(1))
        {
          pt=Point(x/g,y/g);
          return true;
        }
      else
        return false;
    }
  }  
#endif
#endif // CGAL_CARTESIAN_H
#ifdef CGAL_HOMOGENEOUS_H
#ifdef CGAL_LEDA_RATIONAL_H
/*
#ifdef CGAL_LEDA_RATIONAL_H
  bool normalize_coordinates(PointH2<Homogeneous<leda_rational> >& pt)  {
    RT g, x, y, w;
    x = pt.hx();
    y = pt.hy();
    w = pt.hw();
    if (CGAL_NTS is_zero(x) &&  CGAL_NTS is_zero(y)) {
      //g = w;
      pt=Point(x,y,1);
      return true;
    }          
    else {
      g = CGAL_NTS gcd(x, y);
      g = CGAL_NTS gcd(g, w);
      
      pt=Point(x/g,y/g,w/g);
      return true;
    }
  }
#endif
*/
#endif
#ifdef CGAL_LEDA_INTEGER_H
#endif
#endif // CGAL_HOMOGENEOUS_H
#endif // CGAL_POINT_2_H

#ifdef CGAL_ISO_RECTANGLE_2_H
template <class R_> // default does nothing.
inline bool normalize_coordinates(Iso_rectangle_2<R_>& bb)  {
  return false;
}
#ifdef CGAL_CARTESIAN_H
#ifdef CGAL_LEDA_RATIONAL_H
inline bool normalize_coordinates(
              Iso_rectangle_2<Cartesian<leda_rational> >& bb)
{
  //    return normalize(bb[0])&&normalize(bb[2]); 
  // Should be implemented internally in Iso_rectangle, in the near future.
  CGAL::Point_2<Cartesian<leda_rational> > min=bb[0];
  CGAL::Point_2<Cartesian<leda_rational> > max=bb[2];
  if (normalize_coordinates(min)||normalize_coordinates(max)) 
    {
      bb=CGAL::Iso_rectangle_2<Cartesian<leda_rational> >(min,max);
      return true;
    }
  return false;
}
#endif
#ifdef CGAL_LEDA_INTEGER_H
  inline bool normalize_coordinates(
                Iso_rectangle_2<Cartesian<leda_integer> >& bb)
{
    //    return normalize_coordinates(bb[0])&&normalize_coordinates(bb[2]); 
    // Should be implemented internally in Iso_rectangle, in the near future.
    CGAL::Point_2<Cartesian<leda_integer> > min=bb[0];
    CGAL::Point_2<Cartesian<leda_integer> > max=bb[2];
    if (normalize_coordinates(min)||normalize_coordinates(max)) 
      bb=CGAL::Iso_rectangle_2<Cartesian<leda_integer> >(min,max);
    return bb;    
  }
#endif
#endif // CGAL_CARTESIAN_H
#ifdef CGAL_HOMOGENEOUS_H
#ifdef CGAL_LEDA_RATIONAL_H
/*
  inline bool normalize_coordinates(
                Iso_rectangleH2<Homogeneous<leda_rational> >& bb)  {
    return 
    normalize_coordinates((PointH2<Homogeneous<leda_rational> >&) bb[0]) &&
    normalize_coordinates((PointH2<Homogeneous<leda_rational> >&) bb[2]);
  }
*/
#endif
#ifdef CGAL_LEDA_INTEGER_H
#endif
#endif // CGAL_HOMOGENEOUS_H
#endif // CGAL_ISO_RECTANGLE_2_H

CGAL_END_NAMESPACE

#endif // CGAL_PLANAR_MAP_2_NORMALIZE_H










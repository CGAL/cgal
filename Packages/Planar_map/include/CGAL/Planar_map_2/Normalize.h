// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// Author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>

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
    typedef Cartesian<leda_rational> K;
    typedef CGAL::Point_2<K> Point;
    typedef K::RT RT;
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










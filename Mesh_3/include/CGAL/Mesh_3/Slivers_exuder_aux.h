// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent RINEAU

/** This file contains functions from Pierre Alliez, that I have
    "cgalized". */

#ifndef CGAL_MESH_3_SLIVER_EXUDER_AUX_H
#define CGAL_MESH_3_SLIVER_EXUDER_AUX_H

#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Tetrahedron_3.h>

namespace Pierre {

template <typename K>
double
area(const CGAL::Tetrahedron_3<K>& t)
{// This function is (c) 2005 Pierre Alliez
  typedef typename K::Triangle_3 Triangle;
    
  Triangle t1 = Triangle(t[0],t[1],t[2]);
  Triangle t2 = Triangle(t[0],t[1],t[3]);
  Triangle t3 = Triangle(t[2],t[1],t[3]);
  Triangle t4 = Triangle(t[2],t[0],t[3]);
  double a1 = std::sqrt(CGAL_NTS to_double(t1.squared_area()));
  double a2 = std::sqrt(CGAL_NTS to_double(t2.squared_area()));
  double a3 = std::sqrt(CGAL_NTS to_double(t3.squared_area()));
  double a4 = std::sqrt(CGAL_NTS to_double(t4.squared_area()));
  return a1 + a2 + a3 + a4;
}

template <typename K>
double
circumradius(const CGAL::Tetrahedron_3<K>& t)
{ // This function is (c) 2005 Pierre Alliez
  typename K::Point_3 center = circumcenter( t.vertex(0),
                                             t.vertex(1),
                                             t.vertex(2),
                                             t.vertex(3));
  return CGAL::sqrt(CGAL::to_double( squared_distance( center, t.vertex(0))));
}

template <typename K>
double
len(const CGAL::Vector_3<K> &v)
{ // This function is (c) 2005 Pierre Alliez
  return std::sqrt(CGAL::to_double(v*v));
}

double fix_sine(double sine)
{ // This function is (c) 2005 Pierre Alliezx
  if(sine >= 1)
    return 1;
  else
    if(sine <= -1)
      return -1;
    else
      return sine;
}

template <typename K>
double 
angle_rad(const CGAL::Vector_3<K> &u,
          const CGAL::Vector_3<K> &v)
{ // This function is (c) 2005 Pierre Alliez
  // check
  double product = len(u)*len(v);
  if(product == 0.)
    return 0.0;

  // cosine
  double dot = CGAL::to_double(u*v);
  double cosine = dot / product;

  // sine
  typename K::Vector_3 w = CGAL::cross_product(u,v);
  double AbsSine = len(w) / product;

  if(cosine >= 0)
    return std::asin(fix_sine(AbsSine));
  else
    return CGAL_PI-std::asin(fix_sine(AbsSine));
}

template <typename K>
double angle_deg(const CGAL::Vector_3<K> &u,
                 const CGAL::Vector_3<K> &v)
{
  static const double conv = 1.0/CGAL_PI*180.0;

  return conv*angle_rad(u,v);
}

template <typename K>
typename CGAL::Vector_3<K>
normal(const CGAL::Point_3<K>& a,
       const CGAL::Point_3<K>& b,
       const CGAL::Point_3<K>& c)
{
  return CGAL::cross_product(b-a,c-a);
}

// dihedral angle at an edge [vi,vj]
template <typename K>
double dihedral_angle(const CGAL::Tetrahedron_3<K>& t,
                      const int i,
                      const int j)
{
  const CGAL::Vector_3<K> vi = normal(t[(i+1)&3],
                                      t[(i+2)&3],
                                      t[(i+3)&3]);
  const CGAL::Vector_3<K> vj = normal(t[(j+1)&3],
                                      t[(j+2)&3],
                                      t[(j+3)&3]);
  return 180-angle_deg(vi,vj);
}

template <typename K>
double min_dihedral_angle(const CGAL::Tetrahedron_3<K>& t)           
{
  double min = 180;
  
  for(int i = 0; i < 4; ++i)
    for(int j = i + 1; j < 4; j++)
      {
	double a = dihedral_angle(t, i, j);
	if( a < min )
	  min = a;
      }

  std::cerr << min << "\n";
  return min;
}

template <typename K>
double
radius_ratio(const typename CGAL::Tetrahedron_3<K>& t)
{ // This function is (c) 2005 Pierre Alliez
  typename K::FT inradius = 3 * std::abs(t.volume()) / area(t);
  double circ = circumradius(t);
  if(circ == 0)
    return 0;
  else
    return (3 * CGAL::to_double(inradius) / circ);
}

} // end namespace Pierre

#endif // end CGAL_MESH_3_SLIVER_EXUDER_AUX_H

// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef VOLUME_COMPUTATION_H
#define VOLUME_COMPUTATION_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>
#include <cmath>
#include "Element_topo.h"

/*! Compute the area of the triangle p0, p1, p2 (points are 3D points). */
template<typename Point>
typename Point::FT area_of_triangle(const Point& p0,
                                    const Point& p1,
                                    const Point& p2)
{
  typename CGAL::Kernel_traits<Point>::Kernel::Triangle_3 t(p0, p1, p2);
  return std::sqrt(t.squared_area());
}

 /*      3
  *     /|\
  *    0-|-2
  *     \|/
  *      1
  */
template<typename Point>
typename Point::FT signed_volume_of_tetra(const Point& p0, const Point& p1,
                                          const Point& p2, const Point& p3)
{ return CGAL::volume(p0, p1, p2, p3); }

template<typename Point>
typename Point::FT volume_of_tetra(const Point& p0, const Point& p1,
                                   const Point& p2, const Point& p3)
{ return std::abs(signed_volume_of_tetra(p0, p1, p2, p3)); }

/*       4
 *      /|\
 *     0-|-3
 *     | | |
 *     1---2
 */
template<typename Point>
typename Point::FT signed_volume_of_pyramid(const Point& p0, const Point& p1,
                                            const Point& p2, const Point& p3,
                                            const Point& p4)
{ return signed_volume_of_tetra(p0, p1, p2, p4)+
      signed_volume_of_tetra(p0, p2, p3, p4); }

template<typename Point>
typename Point::FT volume_of_pyramid(const Point& p0, const Point& p1,
                                     const Point& p2, const Point& p3,
                                     const Point& p4)
{ return std::abs(signed_volume_of_pyramid(p0, p1, p2, p3, p4)); }

/*      3
 *     /|\
 *    4---5
 *    | | |
 *    | 0 |
 *    |/ \|
 *    1---2
 */
template<typename Point>
typename Point::FT signed_volume_of_prism(const Point& p0, const Point& p1,
                                          const Point& p2, const Point& p3,
                                          const Point& p4, const Point& p5)
{
  return signed_volume_of_tetra(p0, p1, p2, p5) +
      signed_volume_of_tetra(p1, p0, p3, p5) +
      signed_volume_of_tetra(p4, p3, p5, p1); }

template<typename Point>
typename Point::FT volume_of_prism(const Point& p0, const Point& p1,
                                   const Point& p2, const Point& p3,
                                   const Point& p4, const Point& p5)
{ return std::abs(signed_volume_of_prism(p0, p1, p2, p3, p4, p5)); }

/*      7----6
 *     /|   /|
 *    4----5 |
 *    | 3--|-2
 *    |/   |/
 *    0----1
*/
template<typename Point>
typename Point::FT signed_volume_of_hexa(const Point& p0, const Point& p1,
                                  const Point& p2, const Point& p3,
                                  const Point& p4, const Point& p5,
                                  const Point& p6, const Point& p7)
{
  return signed_volume_of_tetra(p0, p1, p3, p4) +
    signed_volume_of_tetra(p3, p1, p2, p6) +
    signed_volume_of_tetra(p6, p5, p4, p1) +
    signed_volume_of_tetra(p6, p4, p7, p3) +
    signed_volume_of_tetra(p1, p3, p4, p6);
}

template<typename Point>
typename Point::FT volume_of_hexa(const Point& p0, const Point& p1,
                                  const Point& p2, const Point& p3,
                                  const Point& p4, const Point& p5,
                                  const Point& p6, const Point& p7)
{ return std::abs(signed_volume_of_hexa(p0, p1, p2, p3, p4, p5, p6, p7)); }

template<typename LCC>
typename LCC::FT signed_volume_of_generic_cell(LCC& lcc,
                                               typename LCC::Dart_handle dh)
{
  typename LCC::FT vol=0;
  typename LCC::Point *p0, *p1, *p2;
  typename LCC::Point p3;
  typename LCC::Point p4=CGAL::ORIGIN; // Used instead barycenter of volume; it is supposed to work whatever the position of this point
  typename LCC::Dart_handle dh2;

  // Iterate through one dart per face of the volume containing dh.
  for (auto it=lcc.template one_dart_per_incident_cell<2,3,2>(dh).begin(),
       itend=lcc.template one_dart_per_incident_cell<2,3,2>(dh).end();
       it!=itend; ++it)
  {
    dh2=lcc.template beta<1,1,1>(it);
    if (dh2==it)
    { // Triangle
      p0=&(lcc.point(it));
      p1=&(lcc.point(lcc.template beta<1>(it)));
      p2=&(lcc.point(lcc.template beta<0>(it)));
      vol+=signed_volume_of_tetra(*p0, *p1, *p2, p4);
    }
    else if  (lcc.template beta<1>(dh2)==it)
    { // Square
      p0=&(lcc.point(it));
      p1=&(lcc.point(lcc.template beta<1>(it)));
      p2=&(lcc.point(lcc.template beta<1,1>(it)));
      vol+=signed_volume_of_tetra(*p0, *p1, *p2, p4);

      p1=p2;
      p2=&(lcc.point(lcc.template beta<0>(it)));
      vol+=signed_volume_of_tetra(*p0, *p1, *p2, p4);
    }
    else
    { // Generic face (>4 edges)
      p3=lcc.template barycenter<2>(it);
      typename LCC::Dart_handle cur=it;
      p0=&(lcc.point(cur));
      do
      {
        cur=lcc.next(cur);
        p1=&(lcc.point(cur));
        vol+=signed_volume_of_tetra(*p0, *p1, p3, p4);
        p0=p1;
      }
      while(cur!=it);
    }
  }
  return vol;
}

template<typename LCC>
typename LCC::FT volume_of_generic_cell(LCC& lcc,
                                        typename LCC::Dart_handle dh)
{ return std::abs(signed_volume_of_generic_cell(lcc, dh)); }

/// @return the volume of the cell
template<typename LCC>
typename LCC::FT signed_volume(LCC& lcc, typename LCC::Dart_handle dh)
{
  typename LCC::Dart_handle sd,d2;
  cell_topo celltopo = Get_cell_topo<LCC, 3>::run(lcc, dh, sd);

  switch(celltopo)
  {
    case TETRAHEDRON:
      return signed_volume_of_tetra(lcc.point(sd),
                                    lcc.point(lcc.template beta<1>(sd)),
                                    lcc.point(lcc.template beta<0>(sd)),
                                    lcc.point(lcc.template beta<2, 0>(sd)));
    case HEXAHEDRON:
      d2=lcc.template beta<2, 1, 1, 2, 1>(sd);
      return signed_volume_of_hexa(lcc.point(sd),
                                   lcc.point(lcc.template beta<1>(sd)),
                                   lcc.point(lcc.template beta<1,1>(sd)),
                                   lcc.point(lcc.template beta<0>(sd)),
                                   lcc.point(d2),
                                   lcc.point(lcc.template beta<0>(d2)),
                                   lcc.point(lcc.template beta<0,0>(d2)),
                                   lcc.point(lcc.template beta<1>(d2)));
    case PRISM:
      d2=lcc.template beta<2, 1, 1, 2>(sd);
      return signed_volume_of_prism(lcc.point(sd),
                                    lcc.point(lcc.template beta<1>(sd)),
                                    lcc.point(lcc.template beta<0>(sd)),
                                    lcc.point(lcc.template beta<1>(d2)),
                                    lcc.point(d2),
                                    lcc.point(lcc.template beta<0>(d2)));
    case PYRAMID:
      return signed_volume_of_pyramid(lcc.point(sd),
                                      lcc.point(lcc.template beta<1>(sd)),
                                      lcc.point(lcc.template beta<1,1>(sd)),
                                      lcc.point(lcc.template beta<0>(sd)),
                                      lcc.point(lcc.template beta<2,0>(sd)));
    case GENERIC_3D:
      return signed_volume_of_generic_cell(lcc, dh);
    default:
      std::cerr<<"Error in signed_volume"<<std::endl;
  }
  return 0;
}

template<typename LCC>
typename LCC::FT volume(LCC& lcc, typename LCC::Dart_handle dh)
{ return std::abs(signed_volume(lcc, dh)); }

#endif // VOLUME_COMPUTATION_H

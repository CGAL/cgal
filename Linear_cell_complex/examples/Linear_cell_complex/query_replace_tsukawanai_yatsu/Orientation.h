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
#ifndef ORIENTATION_H
#define ORIENTATION_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "Element_topo.h"
#include "Volume_computation.h"

/*      3
 *     /|\
 *    0-|-2
 *     \|/
 *      1
 */
template<typename Point>
bool is_positive_tetra(const Point& p0, const Point& p1,
                       const Point& p2, const Point& p3)
{ return CGAL::orientation(p0, p1, p2, p3)==CGAL::POSITIVE; }

/*       4
 *      /|\
 *     0-|-3
 *     | | |
 *     1---2
 */
template<typename Point>
bool is_positive_pyramid(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3,
                         const Point& p4)
{
  CGAL_USE(p3);
  assert(is_positive_tetra(p0, p1, p2, p4)==is_positive_tetra(p0, p2, p3, p4));
  return is_positive_tetra(p0, p1, p2, p4);
}

/*      3
 *     /|\
 *    4---5
 *    | | |
 *    | 0 |
 *    |/ \|
 *    1---2
 */
template<typename Point>
bool is_positive_prism(const Point& p0, const Point& p1,
                       const Point& p2, const Point& p3,
                       const Point& p4, const Point& p5)
{
  CGAL_USE(p3); CGAL_USE(p4);
  assert(is_positive_tetra(p0, p1, p2, p5)==is_positive_tetra(p1, p0, p3, p5));
  assert(is_positive_tetra(p0, p1, p2, p5)==is_positive_tetra(p4, p3, p5, p1));
  return is_positive_tetra(p0, p1, p2, p5);
}

/*      7----6
 *     /|   /|
 *    4----5 |
 *    | 3--|-2
 *    |/   |/
 *    0----1
*/
template<typename Point>
bool is_positive_hexa(const Point& p0, const Point& p1,
                      const Point& p2, const Point& p3,
                      const Point& p4, const Point& p5,
                      const Point& p6, const Point& p7)
{
  CGAL_USE(p2); CGAL_USE(p5); CGAL_USE(p6); CGAL_USE(p7);
  assert(is_positive_tetra(p0, p1, p3, p4)==is_positive_tetra(p3, p1, p2, p6));
  assert(is_positive_tetra(p0, p1, p3, p4)==is_positive_tetra(p6, p5, p4, p1));
  assert(is_positive_tetra(p0, p1, p3, p4)==is_positive_tetra(p6, p4, p7, p3));
  assert(is_positive_tetra(p0, p1, p3, p4)==is_positive_tetra(p1, p3, p4, p6));
  return is_positive_tetra(p0, p1, p3, p4);
}

template<typename LCC>
bool is_positive_generic_cell(LCC& lcc, typename LCC::Dart_handle dh)
{
  // Is is possible to do better? (since virtual tetra are not necessarily correct
  // geometrically...)
  return signed_volume_of_generic_cell(lcc,dh)>0;
}

/// @return test if the 3-cell containing dg has positive orientation
template<typename LCC>
bool is_positive(LCC& lcc, typename LCC::Dart_handle dh)
{
  typename LCC::Dart_handle sd,d2;
  cell_topo celltopo = Get_cell_topo<LCC, 3>::run(lcc, dh, sd);

  switch(celltopo)
  {
    case TETRAHEDRON:
      return is_positive_tetra(lcc.point(sd),
                               lcc.point(lcc.template beta<1>(sd)),
                               lcc.point(lcc.template beta<0>(sd)),
                               lcc.point(lcc.template beta<2, 0>(sd)));
    case HEXAHEDRON:
      d2=lcc.template beta<2, 1, 1, 2, 1>(sd);
      return is_positive_hexa(lcc.point(sd),
                              lcc.point(lcc.template beta<1>(sd)),
                              lcc.point(lcc.template beta<1,1>(sd)),
                              lcc.point(lcc.template beta<0>(sd)),
                              lcc.point(d2),
                              lcc.point(lcc.template beta<0>(d2)),
                              lcc.point(lcc.template beta<0,0>(d2)),
                              lcc.point(lcc.template beta<1>(d2)));
    case PRISM:
      d2=lcc.template beta<2, 1, 1, 2>(sd);
      return is_positive_prism(lcc.point(sd),
                               lcc.point(lcc.template beta<1>(sd)),
                               lcc.point(lcc.template beta<0>(sd)),
                               lcc.point(lcc.template beta<1>(d2)),
                               lcc.point(d2),
                               lcc.point(lcc.template beta<0>(d2)));
    case PYRAMID:
      return is_positive_pyramid(lcc.point(sd),
                                 lcc.point(lcc.template beta<1>(sd)),
                                 lcc.point(lcc.template beta<1,1>(sd)),
                                 lcc.point(lcc.template beta<0>(sd)),
                                 lcc.point(lcc.template beta<2,0>(sd)));
    case GENERIC_3D:
      return is_positive_generic_cell(lcc, dh);
    default:
      std::cerr<<"Error in is_positive"<<std::endl;
  }
  return false;
}

template<typename LCC>
bool check_orientation(LCC& lcc, bool messages=false)
{
  std::size_t nbplus=0, nbminus=0;
  typename LCC::Dart_handle sd;
  bool res=true;

  auto treated_cc=lcc.get_new_mark();
  auto treated_volume=lcc.get_new_mark();
  // std::cout<<"Volumes of each 3-cell:"<<std::endl;
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if(!lcc.is_marked(it, treated_cc))
    {
      std::size_t nbplus_cc=0, nbminus_cc=0;
      for(auto itcc=lcc.template darts_of_cell_basic<4>(it, treated_cc).begin(),
          itccend=lcc.template darts_of_cell_basic<4>(it, treated_cc).end();
          itcc!=itccend; ++itcc)
      {
        lcc.mark(itcc, treated_cc); // Normally not necessary, but safer
        if(!lcc.is_marked(itcc, treated_volume))
        {
          lcc.template mark_cell<3>(itcc, treated_volume);
          if(is_positive(lcc, it)) { ++nbplus_cc; }
          else { ++nbminus_cc; }
        }
      }
      if(nbplus_cc!=0 && nbminus_cc!=0)
      {
        if(messages)
        {
          std::cout<<"[check_orientation]: INCORRECT ORIENTATION: one cc has "
                  <<nbplus_cc<<" positive volumes and "<<nbminus_cc
                  <<" negative ones."<<std::endl;
        }
        res=false;
      }
      nbplus+=nbplus_cc; nbminus+=nbminus_cc;
    }
  }

  if(messages)
  {
    if(nbplus==0)
    {
      if(nbminus==0)
      { std::cout<<"[check_orientation]: empty map."<<std::endl; }
      else
      {
        std::cout<<"[check_orientation]: all the "<<nbminus
                 <<" volumes have negative orientations."<<std::endl;
      }
    }
    else
    {
      if(nbminus==0)
      {
        std::cout<<"[check_orientation]: all the "<<nbplus
                 <<" volumes have positive orientations."<<std::endl;
      }
      else
      {
        std::cout<<"[check_orientation]: INCORRECT ORIENTATION: "
                <<nbplus<<" positive volumes and "<<nbminus
                <<" negative ones."<<std::endl;
      }
    }
  }
  if (nbplus!=0 && nbminus!=0) { res=false; }

  lcc.free_mark(treated_cc);
  lcc.free_mark(treated_volume);
  return res;
}

/// Reorient all the 3-cells to be positive (if true) or negative (if false)
/// @return true iff at least one cell was reorient,
///         false otherwise (the mesh has already the correct orientation)
/// @warning this method does not work if the lcc has no 3-boundary. Indeed,
///          in such a case, the infinite volume of a cc has necessarily the
///          opposite orientation of all the finite volumes in the same cc!!
template<typename LCC>
bool reorient_all(LCC& lcc, bool positive=true)
{
  typename LCC::Dart_handle sd;
  bool res=false;
  bool error=false;

  auto treated_cc=lcc.get_new_mark();
  auto treated_volume=lcc.get_new_mark();
  // std::cout<<"Volumes of each 3-cell:"<<std::endl;
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if(!lcc.is_marked(it, treated_cc))
    {
      std::size_t nbplus_cc=0, nbminus_cc=0;
      for(auto itcc=lcc.template darts_of_cell_basic<4>(it, treated_cc).begin(),
          itccend=lcc.template darts_of_cell_basic<4>(it, treated_cc).end();
          itcc!=itccend; ++itcc)
      {
        lcc.mark(itcc, treated_cc); // Normally not necessary, but safer
        if(!lcc.is_marked(itcc, treated_volume))
        {
          lcc.template mark_cell<3>(itcc, treated_volume);
          if(is_positive(lcc, it)) { ++nbplus_cc; }
          else { ++nbminus_cc; }
        }
      }
      if(nbplus_cc!=0 && nbminus_cc!=0)
      {
        if(!error)
        {
          std::cout<<"[reorient_all]: INCORRECT ORIENTATION: one cc has "
                   <<nbplus_cc<<" positive volumes and "<<nbminus_cc
                   <<" negative ones. It is not possible to directly reorient "
                   <<"the mesh."<<std::endl;
          error=true;
        }
      }
      else if ((nbminus_cc!=0 && positive) || (nbplus_cc!=0 && !positive))
      {
        lcc.reverse_orientation_connected_component(it);
        res=true;
      }
    }
  }

  assert(check_orientation(lcc, false));

  lcc.free_mark(treated_cc);
  lcc.free_mark(treated_volume);
  return res;
}

#endif // ORIENTATION_H

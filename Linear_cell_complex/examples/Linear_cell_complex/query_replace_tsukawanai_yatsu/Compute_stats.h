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
#ifndef COMPUTE_STATS_H
#define COMPUTE_STATS_H

#include "Element_topo.h"
#include "lcc_convexity_test.h"
#include "One_info.h"
#include "Orientation.h"
#include "Volume_computation.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
void display_stats(LCC& lcc, bool with_geom_tests=false)
{
  lcc.display_characteristics(std::cout);
  std::cout<<", valid="<<lcc.is_valid()<<std::endl;

  typename LCC::FT vol=0, sum_vol=0;
  One_info<double> vol_tetra, vol_hexa, vol_prism, vol_pyramid, vol_generic;
  std::size_t nb_0free=0, nb_1free=0, nb_2free=0, nb_3free=0, nb_face_3_free=0;
  std::size_t nbvol=0, nbvolconvex=0;
  cell_topo t;
  typename LCC::Dart_handle sd;

  auto treated=lcc.get_new_mark();
  auto face_treated=lcc.get_new_mark();
  // std::cout<<"Volumes of each 3-cell:"<<std::endl;
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if (lcc.template is_free<0>(it)) { ++nb_0free; }
    if (lcc.template is_free<1>(it)) { ++nb_1free; }
    if (lcc.template is_free<2>(it)) { ++nb_2free; }
    if (lcc.template is_free<3>(it))
    {
      ++nb_3free;
      if (!lcc.is_marked(it, face_treated))
      {
        ++nb_face_3_free;
        lcc.template mark_cell<2>(it, face_treated);
      }
    }
    else
    {
      if (!lcc.is_marked(it, face_treated))
      { lcc.template mark_cell<2>(it, face_treated); }
    }

    if(!lcc.is_marked(it, treated))
    {
      lcc.template mark_cell<3>(it, treated);
      t=Get_cell_topo<LCC,3>::run(lcc, it, sd);
      if(with_geom_tests)
      { vol=signed_volume(lcc, it); }
      switch(t)
      {
        case PRISM:
          //std::cout<<"  PRISM: "<<vol<<std::endl;
          vol_prism.update(vol);
          break;
        case PYRAMID:
          //std::cout<<"  PYRAMID: "<<vol<<std::endl;
          vol_pyramid.update(vol);
          break;
        case TETRAHEDRON:
          //std::cout<<"  TETRAHEDRON: "<<vol<<std::endl;
          vol_tetra.update(vol);
          break;
        case HEXAHEDRON:
          //std::cout<<"  HEXAHEDRON: "<<vol<<std::endl;
          vol_hexa.update(vol);
          break;
        case GENERIC_3D:
          //std::cout<<"  GENERIC_3D: "<<vol<<std::endl;
          vol_generic.update(vol);
          break;
        default:
          //std::cout<<"  ERROR: unknown type: "<<t<<std::endl;
          break;
      }
      if(with_geom_tests)
      {
        if (is_volume_convex(lcc, it)) { ++nbvolconvex; }
        /*std::cout<<"TO DEBUG [display_stats]: vol="<<vol
              <<" generic_vol="<<signed_volume_of_generic_cell(lcc, it)<<std::endl;*/
        sum_vol+=vol;
      }
      ++nbvol;
    }
  }

  std::cout<<"#prisms="<<vol_prism.number_of_elements()
          <<" #pyramids="<<vol_pyramid.number_of_elements()
         <<" #tetrahedra="<<vol_tetra.number_of_elements()
        <<" #hexahedra="<<vol_hexa.number_of_elements()
       <<" #generic="<<vol_generic.number_of_elements()
      <<" (TOTAL="<<nbvol;
  if(with_geom_tests) { std::cout<<", convex="<<nbvolconvex<<")"<<std::endl; }
  else { std::cout<<")"<<std::endl; }

  std::cout<<"#darts free=("<<nb_0free<<", "<<nb_1free<<", "<<nb_2free
          <<", "<<nb_3free<<") #faces 3-free="<<nb_face_3_free<<std::endl;
  if(with_geom_tests)
  {
    std::cout<<"Volume of all the 3-cells: "<<sum_vol<<std::endl;
    std::cout<<"    hexa:"<<vol_hexa<<std::endl;
    std::cout<<"    tetra:"<<vol_tetra<<std::endl;
    std::cout<<"    prism:"<<vol_prism<<std::endl;
    std::cout<<"    pyramid:"<<vol_pyramid<<std::endl;
    std::cout<<"    generic:"<<vol_generic<<std::endl;
  }

  lcc.free_mark(face_treated);
  lcc.free_mark(treated);

  check_orientation(lcc, true);
}
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
void display_volume(LCC& lcc, typename LCC::Dart_handle dh)
{
  typename LCC::Dart_handle sd;
  cell_topo t=Get_cell_topo<LCC,3>::run(lcc, dh, sd);
  switch(t)
  {
    case PRISM:
      std::cout<<"PRISM:";
      break;
    case PYRAMID:
      std::cout<<"PYRAMID:";
      break;
    case TETRAHEDRON:
      std::cout<<"TETRAHEDRON:";
      break;
    case HEXAHEDRON:
      std::cout<<"HEXAHEDRON:";
      break;
    case GENERIC_3D:
      std::cout<<"GENERIC_3D:";
      break;
    default:
      std::cout<<"  ERROR: unknown type: "<<t<<std::endl;
      break;
  }
  std::vector<typename LCC::Point> points;
  for(auto it=lcc.template one_dart_per_incident_cell<0,3>(sd).begin(),
      itend=lcc.template one_dart_per_incident_cell<0,3>(sd).end(); it!=itend; ++it)
  { points.push_back(lcc.point(it)); }
  std::sort(points.begin(), points.end());
  bool first=true;
  for(auto& pt: points)
  {
    if(!first) { std::cout<<"; "; }
    else { std::cout<<"("; first=false; }
    std::cout<<pt;
  }
  std::cout<<")"<<std::endl;
}
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
void display_all_volumes(LCC& lcc)
{
  std::cout<<"*********************************************"<<std::endl;
  for(auto it=lcc.template one_dart_per_cell<3>().begin(),
      itend=lcc.template one_dart_per_cell<3>().end(); it!=itend; ++it)
  { display_volume(lcc, it); }
  std::cout<<"*********************************************"<<std::endl;
}
///////////////////////////////////////////////////////////////////////////////
#endif // COMPUTE_STATS_H

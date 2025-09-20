// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef LCC_READ_FROM_VTK_H
#define LCC_READ_FROM_VTK_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "My_linear_cell_complex_incremental_builder.h"

template<typename LCC>
typename LCC::Dart_handle read_vtk(const std::string& filename, LCC& lcc)
{
  std::ifstream file(filename);
  if(!file.is_open())
  {
    std::cerr<<"[ERROR] read_vtk: cannot open file "<<filename<<std::endl;
    return nullptr;
  }

  typename LCC::Dart_handle res=nullptr;
  My_linear_cell_complex_incremental_builder_3<LCC> ib(lcc);

  std::size_t npoints, ncells, ncells2;
  std::string line,tmp;

  while(line.find("POINTS")==std::string::npos)
  { std::getline(file,line); }
  std::stringstream ss(line);
  std::getline(ss,tmp,' '); // skip POINTS
  ss>>npoints;

  std::vector<typename LCC::Vertex_attribute_handle> points(npoints);
  typename LCC::FT x,y,z;
  for(std::size_t i=0; i<npoints; ++i)
  {
    file>>x>>y>>z;
    points[i]=ib.add_vertex(typename LCC::Point(x, y, z));
  }
  // std::cout<<npoints<<"\n";
  // std::cout<<points[0]<<" "<<points[1]<<" "<<points[2]<<"\n";

  // Read Connectivity
  // read until you find the CELLS  line
  // this is needed because meshes exported from Paraview have sometimes
  // a METADATA section
  while(line.find("CELLS")==std::string::npos)
  { std::getline(file,line); }
  ss=std::stringstream(line);
  std::getline(ss,tmp,' '); // skip CELLS
  ss>>ncells;

  // std::cout<<ncells<<std::endl;

  std::size_t points_per_cell;
  std::vector<std::vector<std::size_t>> faces(ncells);
  for(std::size_t i=0; i<ncells; ++i)
  {
    file>>points_per_cell;
    faces[i].resize(points_per_cell);
    for(std::size_t j=0; j<points_per_cell; ++j)
    { file>>faces[i][j]; }
  }

  while(line.find("CELL_TYPES")==std::string::npos)
  { std::getline(file,line); }
  ss=std::stringstream(line);
  std::getline(ss,tmp,' '); // skip CELL_TYPES
  ss>>ncells2;
  assert(ncells==ncells2);

  std::size_t cell_type;
  for(std::size_t i=0; i<ncells; ++i)
  {
    file>>cell_type;
    //std::cout<<"cell_type "<<cell_type<<"  ";
    switch(cell_type)
    {
      case 10: // TETRA
        res=make_tetrahedron_with_builder(ib,
                                      faces[i][0], faces[i][1],
                                      faces[i][2], faces[i][3]);
        break;
      case 11: // VOXEL (special case in VTK)
        res=make_hexahedron_with_builder(ib, faces[i][0], faces[i][1],
                            faces[i][3], faces[i][2], faces[i][4],
                            faces[i][5], faces[i][7], faces[i][6]);
        break;
      case 12: // HEXA
        res=make_hexahedron_with_builder(ib, faces[i][0], faces[i][1],
                            faces[i][2], faces[i][3], faces[i][4],
                            faces[i][5], faces[i][6], faces[i][7]);
        break;
      case 13: // PRISM (WEDGE in VTK)
        res=make_prism_with_builder(ib, faces[i][0], faces[i][1],
            faces[i][2], faces[i][3], faces[i][4], faces[i][5]);
        break;
      case 14: // PYRAMID
        res=make_pyramid_with_builder(ib, faces[i][0], faces[i][1],
            faces[i][2], faces[i][3], faces[i][4]);
        break;
      case 15: // PENTAGONAL_PRISM
        res=make_pentagonal_prism_with_builder(ib, faces[i][0], faces[i][1],
            faces[i][2], faces[i][3], faces[i][4], faces[i][5],
            faces[i][6], faces[i][7], faces[i][8], faces[i][9]);
        break;
      case 16: // HEXAGONAL_PRISM
        res=make_hexagonal_prism_with_builder(ib, faces[i][0], faces[i][1],
            faces[i][2], faces[i][3], faces[i][4], faces[i][5],
            faces[i][6], faces[i][7], faces[i][8], faces[i][9],
            faces[i][10], faces[i][11]);
        break;
        // TODO: 24 QUADRATIC_TETRA
        //       25 QUADRATIC_HEXAHEDRON
        //       26 QUADRATIC_WEDGE
        //       27 QUADRATIC_PYRAMID

      default:
        std::cerr<<"[ERROR] read_vtk: type "<<cell_type<<" unknown."<<std::endl;
    };
    //std::cout<<std::endl;
  }

  for(auto itv=lcc.vertex_attributes().begin();
      itv!=lcc.vertex_attributes().end(); ++itv)
  {
    if(itv->dart()==nullptr)
    { lcc.erase_vertex_attribute(itv); }
  }
  // Read weights
  // Read until you find the POINTDATA  line
  /*  while(line.find("POINT_DATA") == std::string::npos)
  { std::getline(file,line); }

  std::stringstream ss1(line);

  std::size_t nweights;
  std::getline(ss1,tmp,' '); // skip POINT_DATA
  ss1>>nweights;
  assert(nweights==npoints);
  // skip next two lines
  std::getline(file,line);
  std::getline(file,line);

  std::vector<float> weights(npoints);
  for(std::size_t i = 0 ; i < npoints; i++)
  { file >> weights[i]; }
  std::cout << weights[0] << "\n"; */

  return res;
}

///////////////////////////////////////////////////////////////////////////////
#endif // LCC_READ_FROM_VTK_H

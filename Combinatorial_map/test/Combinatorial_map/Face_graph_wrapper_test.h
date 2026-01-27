// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
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
#ifndef FACE_GRAPH_WRAPPER_TEST_H
#define FACE_GRAPH_WRAPPER_TEST_H

#include <CGAL/Face_graph_wrapper.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

inline bool test_face_graph_wrapper()
{
  bool res=true;

  typedef CGAL::Surface_mesh<CGAL::Simple_cartesian<double>::Point_3> SMesh;
  SMesh m;
  std::ifstream in1(CGAL::data_file_path("meshes/head.off"));
  if (in1.fail())
  {
    std::cout<<"Error: impossible to open 'data/head.off'"<<std::endl;
    return false;
  }
  in1>>m;

  CGAL::Face_graph_wrapper<SMesh> fgw1(m);
  std::vector<unsigned int> cells=fgw1.count_all_cells();
  if (cells[0]!=1487 || cells[1]!=4406 || cells[2]!=2921 ||
      fgw1.number_of_darts()!=8812)
  {
    std::cout<<"Error: incorrect number of cells in test_face_graph_wrapper "
             <<"for Surface_mesh: "
             <<cells[0]<<", "<<cells[1]<<", "<<cells[2]<<", "<<fgw1.number_of_darts()
             <<std::endl;
    res=false;
  }

  typedef CGAL::Polyhedron_3<CGAL::Simple_cartesian<double> > Polyhedron;
  Polyhedron p;
  std::ifstream in2(CGAL::data_file_path("meshes/head.off"));
  if (in2.fail())
  {
    std::cout<<"Error: impossible to open 'data/head.off'"<<std::endl;
    return false;
  }
  in2>>p;
  CGAL::Face_graph_wrapper<Polyhedron> fgw2(p);
  cells=fgw2.count_all_cells();
  if (cells[0]!=1487 || cells[1]!=4406 || cells[2]!=2921 ||
      fgw2.number_of_darts()!=8812)
  {
    std::cout<<"Error: incorrect number of cells in test_face_graph_wrapper "
             <<"for Polyhedron."
             <<cells[0]<<", "<<cells[1]<<", "<<cells[2]<<", "<<fgw2.number_of_darts()
             <<std::endl;
    res=false;
  }

  return res;
}

#endif // FACE_GRAPH_WRAPPER_TEST_H

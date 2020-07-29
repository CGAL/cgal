// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
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
#ifndef GMAP_TEST_INSERTIONS_H
#define GMAP_TEST_INSERTIONS_H

#include <CGAL/Generalized_map_operations.h>

template<typename GMAP>
bool check_number_of_cells_3(GMAP& gmap, unsigned int nbv, unsigned int nbe,
                             unsigned int nbf, unsigned int nbvol,
                             unsigned int nbcc);

template<typename GMAP>
bool test_vertex_insertion(GMAP& gmap)
{
  typename GMAP::Dart_handle d1, d2, d3;

  trace_test_begin();
  d1 = gmap.create_dart();
  gmap.insert_cell_0_in_cell_1(d1);
  if ( !check_number_of_cells_3(gmap, 2, 1, 1, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.create_dart(); gmap.template sew<1>(d1, d1);
  gmap.insert_cell_0_in_cell_1(d1);
  if ( !check_number_of_cells_3(gmap, 2, 1, 1, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_edge();
  gmap.insert_cell_0_in_cell_1(d1);
  if ( !check_number_of_cells_3(gmap, 3, 2, 1, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_edge();
  gmap.template sew<1>(d1, gmap.template alpha<0>(d1));
  gmap.insert_cell_0_in_cell_1(d1);
  if ( !check_number_of_cells_3(gmap, 2, 2, 1, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_edge();
  d2 = gmap.make_edge();
  gmap.template sew<2>(d1, d2);
  gmap.insert_cell_0_in_cell_1(d1);
  if ( !check_number_of_cells_3(gmap, 3, 2, 2, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_edge();
  d2 = gmap.make_edge();
  gmap.template sew<2>(d1, d2);
  gmap.template sew<1>(d1, gmap.template alpha<0>(d1));
  gmap.insert_cell_0_in_cell_1(d1);
  if ( !check_number_of_cells_3(gmap, 2, 2, 2, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_edge();
  d2 = gmap.make_edge();
  gmap.template sew<2>(d1, d2);
  gmap.template sew<1>(d1, gmap.template alpha<0>(d1));
  gmap.template sew<1>(d2, gmap.template alpha<0>(d2));
  gmap.insert_cell_0_in_cell_1(d1);
  if ( !check_number_of_cells_3(gmap, 2, 2, 2, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(3);
  d2 = gmap.template alpha<0,1>(d1);
  d3 = gmap.template alpha<1>(d1);
  gmap.insert_cell_0_in_cell_1(d1);
  gmap.insert_cell_0_in_cell_1(d2);
  gmap.insert_cell_0_in_cell_1(d3);
  if ( !check_number_of_cells_3(gmap, 6, 6, 1, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(3);
  d2 = gmap.make_combinatorial_polygon(3);
  gmap.template sew<3>(d1, d2);
  d2 = gmap.template alpha<0,1>(d1);
  d3 = gmap.template alpha<1>(d1);
  gmap.insert_cell_0_in_cell_1(d1);
  gmap.insert_cell_0_in_cell_1(d2);
  gmap.insert_cell_0_in_cell_1(d3);
  if ( !check_number_of_cells_3(gmap, 6, 6, 1, 2, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(3);
  d2 = gmap.make_combinatorial_polygon(3);
  gmap.template sew<2>(d1, d2);
  gmap.insert_cell_0_in_cell_1(d1);
  if ( !check_number_of_cells_3(gmap, 5, 6, 2, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(3);
  d2 = gmap.make_combinatorial_polygon(3);
  gmap.template sew<2>(d1, d2);
  gmap.insert_cell_0_in_cell_1(gmap.template alpha<0,1>(d1));
  gmap.insert_cell_0_in_cell_1(gmap.template alpha<1>(d1));
  gmap.insert_cell_0_in_cell_1(d1);
  if ( !check_number_of_cells_3(gmap, 7, 8, 2, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_hexahedron();
  d2 = gmap.make_combinatorial_hexahedron();
  gmap.template sew<3>(d1, d2);
  gmap.insert_cell_0_in_cell_1(gmap.template alpha<0,1,0,1>(d1));
  gmap.insert_cell_0_in_cell_1(gmap.template alpha<0,1>(d1));
  gmap.insert_cell_0_in_cell_1(gmap.template alpha<1>(d1));
  gmap.insert_cell_0_in_cell_1(d1);
  if ( !check_number_of_cells_3(gmap, 16, 24, 11, 2, 1) )
    return false;
  gmap.clear();

  return true;
}

template<typename GMAP>
bool test_edge_insertion(GMAP& gmap)
{
  typename GMAP::Dart_handle d1, d2, d3;

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(4);
  gmap.insert_cell_1_in_cell_2(d1, gmap.alpha(d1,0,1,0));
  if ( !check_number_of_cells_3(gmap, 4, 5, 2, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(4);
  d2 = gmap.make_combinatorial_polygon(4);
  gmap.template sew<3>(d1, d2);
  gmap.insert_cell_1_in_cell_2(d1, gmap.alpha(d1,0,1,0));
  if ( !check_number_of_cells_3(gmap, 4, 5, 2, 2, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(4);
  d2 = gmap.make_combinatorial_polygon(4);
  gmap.template sew<2>(d1, d2);
  gmap.insert_cell_1_in_cell_2(d1, gmap.alpha(d1,0,1,0));
  if ( !check_number_of_cells_3(gmap, 6, 8, 3, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.create_dart();
  gmap.insert_dangling_cell_1_in_cell_2(d1);
  if ( !check_number_of_cells_3(gmap, 2, 2, 1, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_edge();
  gmap.template sew<1>(d1, gmap.alpha(d1, 0));
  gmap.insert_dangling_cell_1_in_cell_2(d1);
  if ( !check_number_of_cells_3(gmap, 2, 2, 1, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_edge();
  d2 = gmap.make_edge();
  gmap.template sew<3>(d1, d2);
  gmap.template sew<1>(d1, gmap.alpha(d1, 0));
  gmap.insert_dangling_cell_1_in_cell_2(d1);
  if ( !check_number_of_cells_3(gmap, 2, 2, 1, 2, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(4);
  gmap.insert_dangling_cell_1_in_cell_2(d1);
  if ( !check_number_of_cells_3(gmap, 5, 5, 1, 1, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(4);
  d2 = gmap.make_combinatorial_polygon(4);
  gmap.template sew<3>(d1, d2);
  gmap.insert_dangling_cell_1_in_cell_2(d1);
  if ( !check_number_of_cells_3(gmap, 5, 5, 1, 2, 1) )
    return false;
  gmap.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(4);
  d2 = gmap.make_combinatorial_polygon(4);
  gmap.template sew<2>(d1, d2);
  gmap.insert_dangling_cell_1_in_cell_2(d1);
  if ( !check_number_of_cells_3(gmap, 7, 8, 2, 1, 1) )
    return false;
  gmap.clear();

  return true;
}

template<typename GMAP>
bool test_face_insertion(GMAP& gmap)
{
  typename GMAP::Dart_handle d1, d2, d3;
  std::vector<typename GMAP::Dart_handle> v;

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(4);
  v.push_back(d1);
  v.push_back(gmap.alpha(v[0],0,1));
  v.push_back(gmap.alpha(v[1],0,1));
  v.push_back(gmap.alpha(v[2],0,1));
  gmap.insert_cell_2_in_cell_3(v.begin(),v.end());
  if ( !check_number_of_cells_3(gmap, 4, 4, 2, 1, 1) )
    return false;
  gmap.clear(); v.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_polygon(3);
  d2 = gmap.make_combinatorial_polygon(3);
  gmap.template sew<2>(d1, d2);
  v.push_back(d1);
  v.push_back(gmap.alpha(v[0],0,1));
  v.push_back(gmap.alpha(v[1],0,1));
  gmap.insert_cell_2_in_cell_3(v.begin(),v.end());
  if ( !check_number_of_cells_3(gmap, 4, 5, 3, 2, 1) )
    return false;
  gmap.clear(); v.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_hexahedron();
  d2 = gmap.alpha(d1,2);
  v.push_back(d2);
  v.push_back(gmap.alpha(v[0],0,1,2,1));
  v.push_back(gmap.alpha(v[1],0,1,2,1));
  v.push_back(gmap.alpha(v[2],0,1,2,1));
  gmap.insert_cell_2_in_cell_3(v.begin(),v.end());
  if ( !check_number_of_cells_3(gmap, 8, 12, 7, 2, 1) )
    return false;
  gmap.clear(); v.clear();

  trace_test_begin();
  d1 = gmap.make_combinatorial_hexahedron();
  d2 = gmap.make_combinatorial_hexahedron();
  gmap.template sew<3>(d1,d2);
  d3 = gmap.alpha(d1, 2);
  gmap.template remove_cell<2>(d1);
  v.push_back(d3);
  v.push_back(gmap.alpha(v[0],0,1,2,1));
  v.push_back(gmap.alpha(v[1],0,1,2,1));
  v.push_back(gmap.alpha(v[2],0,1,2,1));
  gmap.insert_cell_2_in_cell_3(v.begin(),v.end());
  if ( !check_number_of_cells_3(gmap, 12, 20, 11, 2, 1) )
    return false;

  GMAP gmap2;
  d1 = gmap2.make_combinatorial_hexahedron();
  d2 = gmap2.make_combinatorial_hexahedron();
  gmap2.template sew<3>(d1,d2);
  if ( !gmap.is_isomorphic_to(gmap2, false, false, false) )
  {
    std::cout<<"Error: gmap and gmap2 are not isomorphic (after insertion/removal).\n";
    assert(false);
    return false;
  }

  if (CGAL::degree<GMAP, 0>(gmap2, d1)!=4)
  {
    std::cout<<"Error: 0-degree is wrong: "<<CGAL::degree<GMAP, 0>(gmap2, d1)<<" instead of 4."<<std::endl;
    assert(false);
    return false;
  }

  if (CGAL::degree<GMAP, 1>(gmap2, d1)!=3)
  {
    std::cout<<"Error: 1-degree is wrong: "<<CGAL::degree<GMAP, 1>(gmap2, d1)<<" instead of 3."<<std::endl;
    assert(false);
    return false;
  }

  if (CGAL::degree<GMAP, 2>(gmap2, d1)!=2)
  {
    std::cout<<"Error: 2-degree is wrong: "<<CGAL::degree<GMAP, 2>(gmap2, d1)<<" instead of 2."<<std::endl;
    assert(false);
    return false;
  }

  if (CGAL::codegree<GMAP, 1>(gmap2, d1)!=2)
  {
    std::cout<<"Error: 1-codegree is wrong: "<<CGAL::codegree<GMAP, 1>(gmap2, d1)<<" instead of 2."<<std::endl;
    assert(false);
    return false;
  }

  if (CGAL::codegree<GMAP, 2>(gmap2, d1)!=4)
  {
    std::cout<<"Error: 2-codegree is wrong: "<<CGAL::codegree<GMAP, 2>(gmap2, d1)<<" instead of 4."<<std::endl;
    assert(false);
    return false;
  }

  if (CGAL::codegree<GMAP, 3>(gmap2, d1)!=6)
  {
    std::cout<<"Error: 3-codegree is wrong: "<<CGAL::codegree<GMAP, 3>(gmap2, d1)<<" instead of 6."<<std::endl;
    assert(false);
    return false;
  }

  gmap.clear(); v.clear();

  return true;
}

#endif // GMAP_TEST_INSERTIONS_H

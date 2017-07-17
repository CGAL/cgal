// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://mcaroli@scm.gforge.inria.fr/svn/cgal/trunk/Periodic_3_triangulation_3/test/Periodic_3_triangulation_3/test_periodic_3_delaunay_3.cpp $
// $Id: test_periodic_3_delaunay_3.cpp 48874 2009-04-23 13:54:38Z mcaroli $
//
//
// Author(s)     : Nico Kruithof
//                 Manuel Caroli

#include <iostream>
#include <fstream>

#include <CGAL/Timer.h>

#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_triangulation_hierarchy_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel          K1;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K1>         PTT1;
typedef CGAL::Periodic_3_triangulation_ds_vertex_base_3<>            DSVB1;
typedef CGAL::Periodic_3_triangulation_ds_cell_base_3<>              DSCB1;
typedef CGAL::Triangulation_vertex_base_3<PTT1,DSVB1>                VBB1;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<VBB1>            VB1;
typedef CGAL::Triangulation_cell_base_3<PTT1,DSCB1>                  CB1;
typedef CGAL::Triangulation_data_structure_3<VB1,CB1>                TDS1;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<PTT1,TDS1>         PDT1;
// Explicit instantiation of the whole class :
template class CGAL::Periodic_3_triangulation_hierarchy_3<PDT1>;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel            K2;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K2>         PTT2;
typedef CGAL::Periodic_3_triangulation_ds_vertex_base_3<>            DSVB2;
typedef CGAL::Periodic_3_triangulation_ds_cell_base_3<>              DSCB2;
typedef CGAL::Triangulation_vertex_base_3<PTT2,DSVB2>                VBB2;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<VBB2>            VB2;
typedef CGAL::Triangulation_cell_base_3<PTT2,DSCB2>                  CB2;
typedef CGAL::Triangulation_data_structure_3<VB2,CB2>                TDS2;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<PTT2,TDS2>         PDT2;
// Explicit instantiation of the whole class :
template class CGAL::Periodic_3_triangulation_hierarchy_3<PDT2>;

#include <CGAL/MP_Float.h>
#include <CGAL/Simple_homogeneous.h>
typedef CGAL::Simple_homogeneous<CGAL::MP_Float>                     K3;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K3>         PTT3;
typedef CGAL::Periodic_3_triangulation_ds_vertex_base_3<>            DSVB3;
typedef CGAL::Periodic_3_triangulation_ds_cell_base_3<>              DSCB3;
typedef CGAL::Triangulation_vertex_base_3<PTT3,DSVB3>                VBB3;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<VBB3>            VB3;
typedef CGAL::Triangulation_cell_base_3<PTT3,DSCB3>                  CB3;
typedef CGAL::Triangulation_data_structure_3<VB3,CB3>                TDS3;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<PTT3,TDS3>         PDT3;
// Explicit instantiation of the whole class :
template class CGAL::Periodic_3_triangulation_hierarchy_3<PDT3>;

#include <CGAL/_test_cls_periodic_3_delaunay_3.h>

int main(int, char**)
{
  CGAL::Timer timer;
  timer.start();
  typedef CGAL::Periodic_3_triangulation_hierarchy_3< PDT1 > P3T3_1;
  _test_cls_periodic_3_delaunay_3( P3T3_1() );

#if ! (defined(_MSC_VER) && defined(_DEBUG))
  typedef CGAL::Periodic_3_triangulation_hierarchy_3< PDT2 > P3T3_2;
  _test_cls_periodic_3_delaunay_3( P3T3_2() );
#endif

  // typedef CGAL::Periodic_3_triangulation_hierarchy_3< PDT3 > P3T3_3;
  // this takes too much time for the test suite.
  //_test_cls_periodic_3_delaunay_3( P3T3_3(), true );

  std::cerr << timer.time() << " sec." << std::endl;
  return 0;
}

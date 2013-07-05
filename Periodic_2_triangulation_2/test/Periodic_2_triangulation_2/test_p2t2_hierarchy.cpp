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
// $URL: svn+ssh://mcaroli@scm.gforge.inria.fr/svn/cgal/trunk/Periodic_2_triangulation_2/test/Periodic_2_triangulation_2/test_periodic_2_delaunay_2.cpp $
// $Id: test_periodic_2_delaunay_2.cpp 48874 2009-04-22 12:54:28Z mcaroli $
//
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>
//                 Manuel Caroli

#include <iostream>
#include <fstream>

#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_triangulation_hierarchy_2.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>

#include "interface_test.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel          K1;
typedef CGAL::Periodic_2_triangulation_traits_2<K1>                  PTT1;
typedef CGAL::Periodic_2_triangulation_vertex_base_2<PTT1>           PVB1;
typedef CGAL::Periodic_2_triangulation_face_base_2<PTT1>             PFB1;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<PVB1>            PHVB1;
typedef CGAL::Triangulation_data_structure_2<PHVB1, PFB1>             Tds1;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<PTT1, Tds1>         PDT1;
// Explicit instantiation of the whole class :
template class CGAL::Periodic_2_triangulation_hierarchy_2<PDT1>;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel            K2;
typedef CGAL::Periodic_2_triangulation_traits_2<K2>                  PTT2;
typedef CGAL::Periodic_2_triangulation_vertex_base_2<PTT2>           DSVB2;
typedef CGAL::Periodic_2_triangulation_face_base_2<PTT2>             DSFB2;
typedef CGAL::Triangulation_vertex_base_2<PTT2, DSVB2>                VBB2;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<VBB2>            VB2;
typedef CGAL::Triangulation_face_base_2<PTT2, DSFB2>                  FB2;
typedef CGAL::Triangulation_data_structure_2<VB2, FB2>                Tds2;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<PTT2, Tds2>         PDT2;
// Explicit instantiation of the whole class :
template class CGAL::Periodic_2_triangulation_hierarchy_2<PDT2>;

#include <CGAL/MP_Float.h>
#include <CGAL/Simple_homogeneous.h>
typedef CGAL::Simple_homogeneous<CGAL::MP_Float>                     K3;
typedef CGAL::Periodic_2_triangulation_traits_2<K3>                  PTT3;
typedef CGAL::Periodic_2_triangulation_vertex_base_2<PTT3>           DSVB3;
typedef CGAL::Periodic_2_triangulation_face_base_2<PTT3>             DSFB3;
typedef CGAL::Triangulation_vertex_base_2<PTT3, DSVB3>                VBB3;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<VBB3>            VB3;
typedef CGAL::Triangulation_face_base_2<PTT3, DSFB3>                  FB3;
typedef CGAL::Triangulation_data_structure_2<VB3, FB3>                Tds3;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<PTT3, Tds3>         PDT3;
// Explicit instantiation of the whole class :
template class CGAL::Periodic_2_triangulation_hierarchy_2<PDT3>;

int main()
{
  std::cout << "New version" << std::endl;
  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  typedef CGAL::Periodic_2_triangulation_hierarchy_2< PDT1 > P2T2_1;
  test<P2T2_1>(false);
  test_delaunay<P2T2_1>();

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  typedef CGAL::Periodic_2_triangulation_hierarchy_2< PDT2 > P2T2_2;
  test<P2T2_2>(true);
  test_delaunay<P2T2_2>();

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  typedef CGAL::Periodic_2_triangulation_hierarchy_2< PDT3 > P2T2_3;
  test<P2T2_3>(true);
  test_delaunay<P2T2_3>();

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  return 0;
}

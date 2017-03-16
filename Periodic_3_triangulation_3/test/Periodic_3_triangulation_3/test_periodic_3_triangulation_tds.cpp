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
// $URL$
// $Id$
//
//
// Author(s)     : Francois Rebufat
//                 Manuel Caroli

#include <CGAL/Timer.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_triangulation_ds_cell_base_3.h>
#include <CGAL/Periodic_3_triangulation_ds_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>

#include <CGAL/_test_cls_periodic_3_tds_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_data_structure_3<
  CGAL::Triangulation_vertex_base_3<
    K,
    CGAL::Periodic_3_triangulation_ds_vertex_base_3<> >,
  CGAL::Triangulation_cell_base_3<
    K,
    CGAL::Periodic_3_triangulation_ds_cell_base_3<> > >       Tds;

// Explicit instantiation :
// template class CGAL::Triangulation_data_structure_3<>;
template class CGAL::Triangulation_data_structure_3<
  CGAL::Triangulation_vertex_base_3<
    K,
    CGAL::Triangulation_ds_vertex_base_3<> >,
  CGAL::Triangulation_cell_base_3<
    K,
    CGAL::Triangulation_ds_cell_base_3<> >
  >;

// just reusing the tests from the T3 package to check whether the
// periodic vertices and cells fulfill the requirements.
int main(int, char**)
{
  CGAL::Timer timer;
  timer.start();
  _test_cls_periodic_3_tds_3(Tds());
  std::cerr << timer.time() << " sec." << std::endl;
  return 0;
}

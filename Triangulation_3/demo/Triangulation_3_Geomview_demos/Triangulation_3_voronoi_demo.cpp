// Copyright (c) 2001, 2002, 2003, 2004, 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

#include <CGAL/basic.h>

#ifndef CGAL_USE_GEOMVIEW
#include <iostream>
int main()
{
  std::cerr << "Geomview doesn't work on this platform,"
               " so this demo doesn't work"
            << std::endl;
  return 0;
}
#else

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>

#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

#include <iostream>

// exact constructions (circumcenter computations) are needed in this
// demo, not only predicates
typedef CGAL::Exact_predicates_exact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_3<K>                 Vb;
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>         TDS;
typedef CGAL::Delaunay_triangulation_3<K, TDS>               Triangulation;
// typedef CGAL::Delaunay_triangulation_3<K> Triangulation;

typedef Triangulation::Point          Point;

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 3, 3, 3));
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.clear();

  Triangulation T;

  std::cout <<"          Inserting points" << std::endl ;
  for (int z=0 ; z<3 ; z++)
    for (int y=0 ; y<3 ; y++)
      for (int x=0 ; x<3 ; x++)
          T.insert(Point(x, y, z));

  T.is_valid(true);

  std::cout <<"          Visualizing T" << std::endl;
  gv.set_wired(true);
  gv << T;

  std::cout <<"          Visualizing the Voronoi edges" << std::endl;
  gv << CGAL::red();
  T.draw_dual(gv);

  char ch;
  std::cout << "Enter any character to quit" << std::endl;
  std::cin >> ch;

  return 0;
}

#endif // CGAL_USE_GEOMVIEW

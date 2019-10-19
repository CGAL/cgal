// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Steve Oudot, Laurent Rineau, Nader Salman

#ifndef CGAL_COMPLEX_2_IN_TRIANGULATION_3_TO_MEDIT_H
#define CGAL_COMPLEX_2_IN_TRIANGULATION_3_TO_MEDIT_H

#include <CGAL/license/Surface_mesher.h>

#include <CGAL/disable_warnings.h>

#include <iomanip>
#include <stack>

namespace CGAL {

template <class C2t3>
void
output_surface_facets_to_medit (std::ostream& os, const C2t3& c2t3)
{
  typedef typename C2t3::Triangulation Tr;
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;

  // Header.
  const Tr& tr = c2t3.triangulation();

  os << "MeshVersionFormatted 1 \n"
     << "Dimension \n"
     << "3\n\n";

  //os << std::setprecision(20);

  // Finite vertices coordinates.
  os << "Vertices\n"
     << tr.number_of_vertices() << " \n";

  std::map<Vertex_handle, int> V;
  int inum = 0;
  for(Finite_vertices_iterator vit = tr.finite_vertices_begin();
      vit != tr.finite_vertices_end();
      ++vit)
  {
    V[vit] = inum++;
    Point p = static_cast<Point>(vit->point());
    os << p.x() << " " << p.y() << " " << p.z() << " 0 \n";
  }


  // Finite facets indices.
  os << "\nTriangles\n"
     << c2t3.number_of_facets() << " \n";

  for( Finite_facets_iterator fit = tr.finite_facets_begin();
                              fit != tr.finite_facets_end(); ++fit)
  {
    if ((*fit).first->is_facet_on_surface((*fit).second)==true)
    {
      for (int i=0; i<4; i++)
        if (i != (*fit).second)
          os << V[(*fit).first->vertex(i)]+1 << " ";

      os << "0 \n"; // without color.
    }
  }
  // Footer
  os << "\nEnd\n";
}

} // end of namespace CGAL

#include <CGAL/enable_warnings.h>

#endif  // CGAL_COMPLEX_2_IN_TRIANGULATION_3_TO_MEDIT_H

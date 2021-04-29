// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// Copyright (c) 2010-2013 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Mael Rouxel-Labbé, Maxime Gimeno, Jane Tournois
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MDS_3_TETRAHEDRON_SOUP_TO_C3T3_H
#define CGAL_MDS_3_TETRAHEDRON_SOUP_TO_C3T3_H

#include <CGAL/license/MDS_3.h>

#include <CGAL/MDS_3/tet_soup_to_c3t3.h>

#include <vector>
#include <array>
#include <map>

namespace CGAL {

  /*!
  \ingroup PkgMDS3Functions

   * @brief 
   * 
   * @tparam TetrahedronRange
   * @tparam Triangulation model of
   * @param tets
   * @param tr
   * 
   * @pre convex volume
  */
  template<typename TetrahedronRange, typename Triangulation>
  void tetrahedron_soup_to_triangulation_3(const TetrahedronRange& tets,
                                           Triangulation& tr)
  {
    typedef Triangulation              Tr;
    typedef typename Tr::Cell_handle   Cell_handle;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Point         Point;

    std::vector<Point> points;
    std::vector<std::array<int, 5> > finite_cells;
    std::map<std::array<int, 3>, typename Tr::Cell::Surface_patch_index> border_facets;
    std::vector<Vertex_handle> vertex_handle_vector;
    std::map<Vertex_handle, int> v2i;

    for (typename TetrahedronRange::value_type tet : tets)
    {
      CGAL_assertion(tet.orientation() != CGAL::NEGATIVE);
      std::array<int, 5> cell;

      Cell_handle hint = Cell_handle();
      for (int i = 0; i < 4; ++i)
      {
        const Point& pi = tet[i];
        typename Tr::Locate_type lt;
        int li, lj;
        hint = tr.locate(pi, lt, li, lj, hint);
        if (lt != Tr::Locate_type::VERTEX)
        {
          points.push_back(pi);
          Vertex_handle newv = tr.insert(pi, lt, hint, li, lj);
          vertex_handle_vector.push_back(newv);

          CGAL_assertion(points.size() == vertex_handle_vector.size());
          v2i.insert(std::make_pair(newv, points.size() - 1));
          cell[i] = static_cast<int>(points.size() - 1);
        }
        else
          cell[i] = v2i.at(hint->vertex(li));
      }
      cell[4] = 1;
      finite_cells.push_back(cell);
    }

    CGAL::MDS_3::build_triangulation(tr, points, finite_cells,
                                     border_facets, vertex_handle_vector);
  }

} //namespace CGAL


#endif // CGAL_MDS_3_TETRAHEDRON_SOUP_TO_C3T3_H

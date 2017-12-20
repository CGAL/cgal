// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//

#ifndef CGAL_TRIANGULATION_3_TO_LCC_H
#define CGAL_TRIANGULATION_3_TO_LCC_H

#include <CGAL/assertions.h>
#include <map>

namespace CGAL {

  /** Convert a given Triangulation_3 into a 3D linear cell complex.
   * @param alcc the used linear cell complex.
   * @param atr the Triangulation_3.
   * @param avol_to_dart a pointer to a std::map associating to each
   *        tetrahedron of atr a corresponding dart in alcc. Not used if NULL.
   * @return A dart incident to the infinite vertex.
   */
  template < class LCC, class Triangulation >
  typename LCC::Dart_handle import_from_triangulation_3
  (LCC& alcc, const Triangulation &atr,
   std::map<typename Triangulation::Cell_handle,
            typename LCC::Dart_handle >* avol_to_dart=NULL)
  {
    CGAL_static_assertion( LCC::dimension>=3 && LCC::ambient_dimension==3 );

    // Case of empty triangulations.
    if (atr.number_of_vertices() == 0) return LCC::null_handle;

    // Check the dimension.
    if (atr.dimension() != 3) return LCC::null_handle;
    CGAL_assertion(atr.is_valid());

    typedef typename Triangulation::Vertex_handle    TVertex_handle;
    typedef typename Triangulation::Vertex_iterator  TVertex_iterator;
    typedef typename Triangulation::Cell_iterator    TCell_iterator;
    typedef typename std::map
      < TCell_iterator, typename LCC::Dart_handle >::iterator itmap_tcell;

    // Create vertices in the map and associate in a map
    // TVertex_handle and vertices in the map.
    std::map< TVertex_handle, typename LCC::Vertex_attribute_handle > TV;
    for (TVertex_iterator itv = atr.vertices_begin();
         itv != atr.vertices_end(); ++itv)
    {
      TV[itv] = alcc.create_vertex_attribute(itv->point());
    }

    // Create the tetrahedron and create a map to link Cell_iterator
    // and tetrahedron.
    TCell_iterator it;

    std::map<typename Triangulation::Cell_handle, typename LCC::Dart_handle> TC;
    std::map<typename Triangulation::Cell_handle, typename LCC::Dart_handle>*
      mytc = (avol_to_dart==NULL?&TC:avol_to_dart);

    itmap_tcell maptcell_it;

    typename LCC::Dart_handle res=LCC::null_handle, dart=LCC::null_handle;
    typename LCC::Dart_handle cur=LCC::null_handle, neighbor=LCC::null_handle;

    for (it = atr.cells_begin(); it != atr.cells_end(); ++it)
    {
      /*     if (it->vertex(0) != atr.infinite_vertex() &&
             it->vertex(1) != atr.infinite_vertex() &&
             it->vertex(2) != atr.infinite_vertex() &&
             it->vertex(3) != atr.infinite_vertex())
      */
      {
        res = alcc.make_tetrahedron(TV[it->vertex(0)],
                                    TV[it->vertex(1)],
                                    TV[it->vertex(2)],
                                    TV[it->vertex(3)]);

        if ( dart==LCC::null_handle )
        {
          if ( it->vertex(0) == atr.infinite_vertex() )
            dart = res;
          else if ( it->vertex(1) == atr.infinite_vertex() )
            dart = alcc.next(res);
          else if ( it->vertex(2) == atr.infinite_vertex() )
            dart = alcc.previous(res);
          else if ( it->vertex(3) == atr.infinite_vertex() )
            dart = alcc.previous(alcc.template opposite<2>(res));
        }

        for (unsigned int i = 0; i < 4; ++i)
        {
          switch (i)
          {
          case 0: cur = alcc.template opposite<2>(alcc.next(res)); break;
          case 1: cur = alcc.template opposite<2>(alcc.previous(res)); break;
          case 2: cur = alcc.template opposite<2>(res); break;
          case 3: cur = res; break;
          }

          maptcell_it = mytc->find(it->neighbor(i));
          if (maptcell_it != mytc->end())
          {
            switch (atr.mirror_index(it,i) )
            {
            case 0: neighbor = alcc.template opposite<2>(alcc.next(maptcell_it->second));
              break;
            case 1: neighbor = alcc.template opposite<2>(alcc.previous(maptcell_it->second));
              break;
            case 2: neighbor = alcc.template opposite<2>(maptcell_it->second); break;
            case 3: neighbor = maptcell_it->second; break;
            }
            while (alcc.vertex_attribute(neighbor) !=
                   alcc.vertex_attribute(alcc.other_extremity(cur)) )
              neighbor = alcc.next(neighbor);
            alcc.template topo_sew<3>(cur, alcc.other_orientation(neighbor));
          }
        }
        (*mytc)[it] = res;
      }
    }
    CGAL_assertion(dart!=LCC::null_handle);
    return dart;
  }

} // namespace CGAL

#endif // CGAL_TRIANGULATION_3_TO_LCC_H

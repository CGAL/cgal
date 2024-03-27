// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//

#ifndef CGAL_TRIANGULATION_2_TO_LCC_H
#define CGAL_TRIANGULATION_2_TO_LCC_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/assertions.h>
#include <map>

namespace CGAL {

  /** Convert a given Triangulation_2 into a 2D linear cell complex.
   * @param alcc the used linear cell complex.
   * @param atr the Triangulation_2.
   * @param aface_to_dart a pointer to a std::map associating to each
   *        triangle of atr a corresponding dart in alcc. Not used if nullptr.
   * @return A dart incident to the infinite vertex.
   */
  template < class LCC, class Triangulation >
  typename LCC::Dart_descriptor import_from_triangulation_2
  (LCC& alcc, const Triangulation &atr,
   std::map<typename Triangulation::Face_handle,
            typename LCC::Dart_descriptor >* aface_to_dart=nullptr)
  {
    static_assert( LCC::dimension>=2 && LCC::ambient_dimension==2 );

    // Case of empty triangulations.
    if (atr.number_of_vertices()==0) return LCC::null_descriptor;

    // Check the dimension.
    if (atr.dimension()!=2) return LCC::null_descriptor;
    CGAL_assertion(atr.is_valid());

    typedef typename Triangulation::Vertex_handle         TVertex_handle;
    typedef typename Triangulation::All_vertices_iterator TVertex_iterator;
    typedef typename Triangulation::All_faces_iterator    TFace_iterator;
    typedef typename std::map
      < TFace_iterator, typename LCC::Dart_descriptor >::iterator itmap_tcell;

    // Create vertices in the map and associate in a map
    // TVertex_handle and vertices in the map.
    std::map< TVertex_handle, typename LCC::Vertex_attribute_descriptor > TV;
    for (TVertex_iterator itv = atr.all_vertices_begin();
         itv != atr.all_vertices_end(); ++itv)
    {
      TV[itv] = alcc.create_vertex_attribute(itv->point());
    }

    // Create the triangles and create a map to link Cell_iterator
    // and triangles.
    TFace_iterator it;

    std::map<typename Triangulation::Face_handle, typename LCC::Dart_descriptor> TC;
    std::map<typename Triangulation::Face_handle, typename LCC::Dart_descriptor>*
      mytc = (aface_to_dart==nullptr?&TC:aface_to_dart);

    itmap_tcell maptcell_it;

    typename LCC::Dart_descriptor res=LCC::null_descriptor, dart=LCC::null_descriptor;
    typename LCC::Dart_descriptor cur=LCC::null_descriptor, neighbor=LCC::null_descriptor;

    for (it = atr.all_faces_begin(); it != atr.all_faces_end(); ++it)
    {
      /*     if (it->vertex(0) != atr.infinite_vertex() &&
             it->vertex(1) != atr.infinite_vertex() &&
             it->vertex(2) != atr.infinite_vertex() &&
             it->vertex(3) != atr.infinite_vertex())
      */
      {
        res = alcc.make_triangle(TV[it->vertex(0)],
                                 TV[it->vertex(1)],
                                 TV[it->vertex(2)]);

        if ( dart==LCC::null_descriptor )
        {
          if ( it->vertex(0) == atr.infinite_vertex() )
            dart = res;
          else if ( it->vertex(1) == atr.infinite_vertex() )
            dart = alcc.next(res);
          else if ( it->vertex(2) == atr.infinite_vertex() )
            dart = alcc.previous(res);
        }

        for (unsigned int i=0; i<3; ++i)
        {
          switch (i)
          {
          case 0: cur = alcc.next(res); break;
          case 1: cur = alcc.previous(res); break;
          case 2: cur = res; break;
          }

          maptcell_it = mytc->find(it->neighbor(i));
          if (maptcell_it != mytc->end())
          {
            switch (atr.mirror_index(it,i) )
            {
            case 0: neighbor = alcc.next(maptcell_it->second);
              break;
            case 1: neighbor = alcc.previous(maptcell_it->second);
              break;
            case 2: neighbor = maptcell_it->second; break;
            }
            alcc.template set_opposite<2>(cur, neighbor);
          }
        }
        (*mytc)[it] = res;
      }
    }

    CGAL_assertion(dart!=LCC::null_descriptor);
    return dart;
  }

} // namespace CGAL

#endif // CGAL_TRIANGULATION_2_TO_LCC_H

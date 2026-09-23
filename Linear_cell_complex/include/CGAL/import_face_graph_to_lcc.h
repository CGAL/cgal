// Copyright (c) 2011-2026 CNRS and LIRIS' Establishments (France).
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

#ifndef CGAL_IMPORT_FACE_GRAPH_TO_LCC_H
#define CGAL_IMPORT_FACE_GRAPH_TO_LCC_H

//~ #include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/iterator.h>


namespace CGAL {


/*!
\ingroup PkgLinearCellComplexConstructions

imports `pmesh` into `lcc`. The polygon mesh is added in `lcc`, existing darts are not modified. Returns a dart created during the import.
\pre \link GenericMap::dimension `LCC::dimension`\endlink \f$ \geq\f$ 2 and \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3.

\tparam LCC a model of `LinearCellComplex`
\tparam PolygonMesh a model of `FaceListGraph`

\sa `CGAL::read_plane_graph_in_lcc<LCC>`
\sa `CGAL::triangulation_3_to_lcc<LCC,Triangulation>`
*/
template<class PolygonMesh, class LCC>
typename LCC::Dart_descriptor import_face_graph_to_lcc(const PolygonMesh &pmesh,
                                                       LCC& lcc)
{
  static_assert( LCC::dimension>=2 && LCC::ambient_dimension==3 );

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  Halfedge_handle;
  typedef typename boost::graph_traits<PolygonMesh>::face_iterator   Facet_iterator;
  typedef Halfedge_around_face_circulator<PolygonMesh> HF_circulator;

  typedef std::map < Halfedge_handle, typename LCC::Dart_descriptor>
    Halfedge_handle_map;
  typedef typename Halfedge_handle_map::iterator itmap_hds;
  Halfedge_handle_map TC;

  itmap_hds it;
  typename LCC::Dart_descriptor d = LCC::null_descriptor, prev = LCC::null_descriptor;
  typename LCC::Dart_descriptor firstFacet = LCC::null_descriptor, firstAll = LCC::null_descriptor;

  // First traversal to build the darts and link them.
  for (Facet_iterator i = faces(pmesh).begin(); i != faces(pmesh).end(); ++i)
  {
    HF_circulator j(halfedge(*i,pmesh), pmesh), done(j);
    prev = LCC::null_descriptor;
    do
    {
      d = lcc.make_half_edge();
      TC[*j] = d;

      if (prev != LCC::null_descriptor) lcc.set_next(prev, d);
      else firstFacet = d;

      if (!is_border(opposite(*j, pmesh), pmesh))
      {
        it = TC.find(opposite(*j,pmesh));
        if (it != TC.end())
          lcc.template set_opposite<2>(d, it->second);
      }
      prev = d;
    }
    while (++j != done);
    lcc.set_next(prev, firstFacet);
    if (firstAll == LCC::null_descriptor) firstAll = firstFacet;
  }

  // Second traversal to update the geometry.
  // We run one again through the facets of the HDS.
  auto vpm = get(CGAL::vertex_point, pmesh);
  for (Facet_iterator i = faces(pmesh).begin(); i != faces(pmesh).end(); ++i)
  {
    HF_circulator j(halfedge(*i,pmesh), pmesh), done(j);
    do
    {
      d = TC[*j]; // Get the dart associated to the Halfedge
      if (lcc.vertex_attribute(d)==LCC::null_descriptor)
      {
        lcc.set_vertex_attribute
          (d, lcc.create_vertex_attribute(get(vpm, target(opposite(*j, pmesh), pmesh))));
      }
    }
    while (++j != done);
  }

  return firstAll;
}

#ifdef DOXYGEN_RUNNING
/*!
\ingroup PkgLinearCellComplexConstructions

@deprecated This function is deprecated. Users should instead use `CGAL::import_face_graph_to_lcc()`


Imports `apoly`into `lcc`. Objects are added in `lcc`, existing darts are not modified. Returns a dart created during the import.
\pre \link GenericMap::dimension `LCC::dimension`\endlink \f$ \geq\f$ 2 and \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3.

\tparam LCC a model of `LinearCellComplex`
\tparam PolygonMesh a model of `FaceGraph`

\sa `CGAL::read_plane_graph_in_lcc<LCC>`
\sa `CGAL::triangulation_3_to_lcc<LCC,Triangulation>`
*/
template<class LCC,class PolygonMesh>
typename LCC::Dart_descriptor import_from_polyhedron_3(LCC& lcc, const PolygonMesh &apoly);


/*!
\ingroup PkgLinearCellComplexConstructions

@deprecated This function is deprecated. Users should instead use `CGAL::import_face_graph_to_lcc()`


Imports `apoly`into `lcc`. Objects are added in `lcc`, existing darts are not modified. Returns a dart created during the import.
\pre \link GenericMap::dimension `LCC::dimension`\endlink \f$ \geq\f$ 2 and \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3.

\tparam LCC a model of `LinearCellComplex`
\tparam PolygonMesh a model of `FaceGraph`

\sa `CGAL::read_plane_graph_in_lcc<LCC>`
\sa `CGAL::triangulation_3_to_lcc<LCC,Triangulation>`
*/
template<class LCC,class PolygonMesh>
typename LCC::Dart_descriptor polyhedron_3_to_lcc(LCC& lcc, const PolygonMesh &apoly);

#endif

} // end of CGAL namespace

#endif

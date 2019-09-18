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

#ifndef CGAL_POLYHEDRON_3_TO_LCC_H
#define CGAL_POLYHEDRON_3_TO_LCC_H

#include <CGAL/assertions.h>
#include <iostream>
#include <map>
#include <CGAL/Polyhedron_3.h>

namespace CGAL {

  /** Import a given Polyhedron_3 into a Linear_cell_complex.
   * @param alcc the linear cell complex where Polyhedron_3 will be converted.
   * @param apoly the Polyhedron.
   * @return A dart created during the convertion.
   */
  template< class LCC, class Polyhedron >
  typename LCC::Dart_handle import_from_polyhedron_3(LCC& alcc,
                                                     const Polyhedron &apoly)
  {
    CGAL_static_assertion( LCC::dimension>=2 && LCC::ambient_dimension==3 );

    typedef typename Polyhedron::Halfedge_const_handle  Halfedge_handle;
    typedef typename Polyhedron::Facet_const_iterator   Facet_iterator;
    typedef typename Polyhedron::Halfedge_around_facet_const_circulator
      HF_circulator;

    typedef std::map < Halfedge_handle, typename LCC::Dart_handle>
      Halfedge_handle_map;
    typedef typename Halfedge_handle_map::iterator itmap_hds;
    Halfedge_handle_map TC;

    itmap_hds it;
    typename LCC::Dart_handle d = LCC::null_handle, prev = LCC::null_handle;
    typename LCC::Dart_handle firstFacet = LCC::null_handle, firstAll = LCC::null_handle;

    // First traversal to build the darts and link them.
    for (Facet_iterator i = apoly.facets_begin(); i != apoly.facets_end(); ++i)
    {
      HF_circulator j = i->facet_begin();
      prev = LCC::null_handle;
      do
      {
        d = alcc.make_half_edge();
        TC[j] = d;

        if (prev != LCC::null_handle) alcc.set_next(prev, d);
        else firstFacet = d;

        if (!j->opposite()->is_border())
        {
          it = TC.find(j->opposite());
          if (it != TC.end())
            alcc.template set_opposite<2>(d, it->second);
        }
        prev = d;
      }
      while (++j != i->facet_begin());
      alcc.set_next(prev, firstFacet);
      if (firstAll == LCC::null_handle) firstAll = firstFacet;
    }

    // Second traversal to update the geometry.
    // We run one again through the facets of the HDS.
    for (Facet_iterator i = apoly.facets_begin(); i != apoly.facets_end(); ++i)
    {
      HF_circulator j = i->facet_begin();
      do
      {
        d = TC[j]; // Get the dart associated to the Halfedge
        if (alcc.vertex_attribute(d)==LCC::null_handle)
        {
          alcc.set_vertex_attribute
            (d, alcc.create_vertex_attribute(j->opposite()->vertex()->point()));
        }
      }
      while (++j != i->facet_begin());
    }
    return firstAll;
  }

  /** Convert a Polyhedron_3 read into a flux into 3D linear cell complex.
   * @param alcc the linear cell complex where Polyhedron_3 will be converted.
   * @param ais the istream where read the Polyhedron_3.
   * @return A dart created during the convertion.
   */
  template < class LCC >
  typename LCC::Dart_handle
  import_from_polyhedron_3_flux(LCC& alcc, std::istream& ais)
  {
    if (!ais.good())
    {
      std::cout << "Error reading flux." << std::endl;
      return LCC::null_handle;
    }
    CGAL::Polyhedron_3<typename LCC::Traits> P;
    ais >> P;
    return import_from_polyhedron_3<LCC, CGAL::Polyhedron_3
                                    <typename LCC::Traits> > (alcc, P);
  }

} // namespace CGAL

#endif // CGAL_IMPORT_FROM_POLYHEDRON_3_H

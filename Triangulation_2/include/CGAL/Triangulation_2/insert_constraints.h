// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Andreas Fabri, Mariette Yvinec

#ifndef CGAL_INTERNAL_TRIANGULATION_2_IMSERT_CONSTRAINTS_H
#define CGAL_INTERNAL_TRIANGULATION_2_IMSERT_CONSTRAINTS_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <vector>
#include <iterator>

namespace CGAL {
  namespace internal {



    template <class T, class IndicesIterator>
    std::size_t insert_constraints( T& t,
                                    const std::vector<typename T::Point>& points,
                                    IndicesIterator indices_first,
                                    IndicesIterator indices_beyond )
  {
    if(indices_first == indices_beyond){
      return 0;
    }
    typedef typename T::Vertex_handle Vertex_handle;
    typedef typename T::Face_handle Face_handle;
    typedef typename T::Geom_traits Geom_traits;
    typedef typename T::Point Point;
    typedef std::vector<std::size_t> Vertex_indices;
    typedef std::vector<Vertex_handle> Vertices;

    Vertex_indices vertex_indices;
    vertex_indices.resize(points.size());

    std::copy(boost::counting_iterator<std::ptrdiff_t>(0),
              boost::counting_iterator<std::ptrdiff_t>(points.size()),
              std::back_inserter(vertex_indices));

    typename T::size_type n = t.number_of_vertices();
    CGAL::Spatial_sort_traits_adapter_2<
      Geom_traits,
      typename Pointer_property_map<Point>::const_type >
        sort_traits(make_property_map(points),t.geom_traits());

    spatial_sort(vertex_indices.begin(), vertex_indices.end(), sort_traits);

    Vertices vertices;
    vertices.resize(points.size());

    Face_handle hint;
    for(typename Vertex_indices::const_iterator
          it_pti = vertex_indices.begin(), end = vertex_indices.end();
          it_pti != end; ++it_pti)
    {
      vertices[*it_pti] = t.insert(points[*it_pti], hint);
      hint=vertices[*it_pti]->face();
    }

    for(IndicesIterator it_cst=indices_first, end=indices_beyond;
        it_cst!=end; ++it_cst)
    {
      Vertex_handle v1 = vertices[it_cst->first];
      Vertex_handle v2 = vertices[it_cst->second];
      if(v1 != v2) t.insert_constraint(v1, v2);
    }

    return t.number_of_vertices() - n;
  }



    template <class T,class ConstraintIterator>
    std::size_t insert_constraints(T& t,
                                   ConstraintIterator first,
                                   ConstraintIterator beyond)
  {
    typedef typename T::Point Point;
    typedef typename T::Point Point;
    std::vector<Point> points;
    for (ConstraintIterator s_it=first; s_it!=beyond; ++s_it)
    {
      points.push_back( T::get_source(*s_it) );
      points.push_back( T::get_target(*s_it) );
    }

    std::vector< std::pair<std::size_t, std::size_t> > segment_indices;
    std::size_t nb_segments = points.size() / 2;
    segment_indices.reserve( nb_segments );
    for (std::size_t k=0; k < nb_segments; ++k)
      segment_indices.push_back( std::make_pair(2*k,2*k+1) );

    return insert_constraints( t,
                               points,
                               segment_indices.begin(),
                               segment_indices.end() );
  }


  }
}

#endif // CGAL_INTERNAL_TRIANGULATION_2_IMSERT_CONSTRAINTS_H

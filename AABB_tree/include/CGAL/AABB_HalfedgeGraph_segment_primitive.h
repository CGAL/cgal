// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sebastien Loriot
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_AABB_HALFEDGEGRAPH_TRIANGLE_PRIMITIVE_H
#define CGAL_AABB_HALFEDGEGRAPH_TRIANGLE_PRIMITIVE_H

#include <CGAL/AABB_primitive.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_3_property_map.h>

#include <iterator>
#include <boost/mpl/and.hpp>
#include <CGAL/is_iterator.h>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

namespace CGAL {
  
template < class HalfedgeGraph,
           class cache_datum=Tag_false,
           class Id_=typename boost::graph_traits<HalfedgeGraph>::edge_descriptor
           >
class AABB_HalfedgeGraph_segment_primitive : public AABB_primitive< Id_,
                                                                    Segment_from_edge_descriptor_property_map<HalfedgeGraph>,
                                                                    Source_point_from_edge_descriptor<HalfedgeGraph>,
                                                                    Tag_true,
                                                                    cache_datum >
{
  typedef Segment_from_edge_descriptor_property_map<HalfedgeGraph>  Triangle_property_map;
  typedef Source_point_from_edge_descriptor<HalfedgeGraph> Point_property_map;
  
  typedef AABB_primitive< Id_,
                          Triangle_property_map,
                          Point_property_map,
                          Tag_true,
                          cache_datum > Base;
  
public:
  // constructors
  template <class Iterator>
  AABB_HalfedgeGraph_segment_primitive(Iterator it, const HalfedgeGraph& graph)
    : Base( Id_(*it),
            Triangle_property_map(&graph),
            Point_property_map(&graph) ){}

  //for backward-compatibility with AABB_polyhedron_segment_primitive
  AABB_HalfedgeGraph_segment_primitive(Id_ id)
    : Base( id,
            Triangle_property_map(NULL),
            Point_property_map(NULL) ){}
              
  static typename Base::Extra_data construct_primitive_data( const HalfedgeGraph& graph )
  {
    return Base::construct_primitive_data(Triangle_property_map(&graph), Point_property_map(&graph));
  }
  
  //for backward-compatibility with AABB_polyhedron_segment_primitive
  static typename Base::Extra_data construct_primitive_data()
  {
    return Base::construct_primitive_data(Triangle_property_map(NULL), Point_property_map(NULL));
  }
};

}  // end namespace CGAL


#endif // CGAL_AABB_HALFEDGEGRAPH_TRIANGLE_PRIMITIVE_H


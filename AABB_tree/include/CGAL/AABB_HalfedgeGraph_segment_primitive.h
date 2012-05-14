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
  
namespace internal{
  template <class T,bool is_iterator_=is_iterator<T>::value>
  struct Extract_value_type_of_iterator{
    struct value_type{};
  };

  template <class T>
  struct Extract_value_type_of_iterator<T,true>{
    typedef typename std::iterator_traits<T>::value_type value_type;
  };
}
  
template < class HalfedgeGraph,
           class cache_datum=Tag_false,
           class Id_=typename boost::graph_traits<HalfedgeGraph>::edge_descriptor
           >
class AABB_HalfedgeGraph_segment_primitive : public AABB_primitive< Id_,
                                                                    Segment_from_edge_descriptor_property_map<HalfedgeGraph>,
                                                                    Source_point_from_edge_descriptor<HalfedgeGraph>,
                                                                    Tag_false,
                                                                    cache_datum >
{
  typedef Segment_from_edge_descriptor_property_map<HalfedgeGraph>  Triangle_property_map;
  typedef Source_point_from_edge_descriptor<HalfedgeGraph> Point_property_map;
  
  typedef AABB_primitive< Id_,
                          Triangle_property_map,
                          Point_property_map,
                          Tag_false,
                          cache_datum > Base;
  
public:
  // constructors
  AABB_HalfedgeGraph_segment_primitive(Id_ it) : Base(it){}
  //the enable_if is required here so that the first overload is chosen when an Edge_iterator is provided
  template <class Iterator>
  AABB_HalfedgeGraph_segment_primitive(Iterator it,
                                       typename boost::enable_if<
                                          boost::is_convertible<
                                            typename internal::Extract_value_type_of_iterator<Iterator>::value_type,
                                            std::pair<HalfedgeGraph*,Id_>
                                          >
                                        >::type* = NULL )
    : Base( it->second,
            Triangle_property_map((it->first)),
            Point_property_map((it->first)) ){}
};

}  // end namespace CGAL


#endif // CGAL_AABB_HALFEDGEGRAPH_TRIANGLE_PRIMITIVE_H


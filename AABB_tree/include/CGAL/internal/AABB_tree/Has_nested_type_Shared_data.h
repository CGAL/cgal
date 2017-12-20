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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_INTERNAL_AABB_TREE_HAS_NESTED_TYPE_SHARED_DATA_H
#define CGAL_INTERNAL_AABB_TREE_HAS_NESTED_TYPE_SHARED_DATA_H

#include <CGAL/license/AABB_tree.h>


#include <boost/mpl/has_xxx.hpp>
#include <CGAL/tags.h>

namespace CGAL{

namespace internal{

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Shared_data,Shared_data,false)

// Utility class used by AABB_face_graph_triangle_primitive and AABB_halfedge_graph_segment_primitive
// to implement the Construct_shared_data static function.
template<class Graph, class Base, class ObjectPropertyMap, class PointPropertyMap, class HasSharedDataTag>
struct Cstr_shared_data;

template<class Graph, class Base, class ObjectPropertyMap, class PointPropertyMap>
struct Cstr_shared_data<Graph, Base, ObjectPropertyMap, PointPropertyMap, ::CGAL::Tag_true>
{
  typedef typename Base::Shared_data Shared_data;
  static Shared_data construct_shared_data(Graph& graph)
  {
    return Base::construct_shared_data(ObjectPropertyMap(&graph), PointPropertyMap(&graph));
  }

  template <class VertexPmap>
  static Shared_data construct_shared_data(Graph& graph, const VertexPmap& vpm)
  {
    return Base::construct_shared_data(ObjectPropertyMap(&graph, vpm), PointPropertyMap(&graph, vpm));
  }
};

template<class Graph, class Base, class ObjectPropertyMap, class PointPropertyMap>
struct Cstr_shared_data<Graph, Base, ObjectPropertyMap, PointPropertyMap, ::CGAL::Tag_false>
{
  typedef void* Shared_data;
  static Shared_data construct_shared_data(Graph&)
  {
    return NULL;
  }

  template <class VertexPmap>
  static Shared_data construct_shared_data(Graph&, VertexPmap)
  {
    return NULL;
  }
};

} } //namespace CGAL

#endif //CGAL_INTERNAL_AABB_TREE_HAS_NESTED_TYPE_SHARED_DATA_H

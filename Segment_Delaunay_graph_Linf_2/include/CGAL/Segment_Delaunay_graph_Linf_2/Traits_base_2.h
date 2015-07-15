// Copyright (c) 2015  Università della Svizzera italiana.
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
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_TRAITS_BASE_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_TRAITS_BASE_2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2/basic.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Constructions_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Traits_base_2.h>

namespace CGAL {

//-----------------------------------------------------------------------
// the Traits class
//-----------------------------------------------------------------------

template<class R, class MTag, class ITag>
class Segment_Delaunay_graph_Linf_traits_base_2
 : public Segment_Delaunay_graph_traits_base_2<R, MTag, ITag>
{
public:
  //-----------------------------------------------------------------------
  //                  TYPE DEFINITIONS
  //-----------------------------------------------------------------------

  // BASIC TYPES
  //------------
  
  typedef Segment_Delaunay_graph_traits_base_2<R, MTag, ITag> Base;

  typedef typename Base::Kernel Kernel;
  typedef Kernel   K;

  typedef typename Kernel::Iso_rectangle_2        Iso_rectangle_2;
  typedef typename Kernel::Direction_2            Direction_2;
  typedef typename Kernel::Vector_2               Vector_2;

  typedef typename Kernel::Sign                   Sign;
  typedef typename Kernel::Boolean                Boolean;

  // CONSTRUCTIONS
  //--------------
  // vertex and bisector
  typedef
  CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::Construct_svd_vertex_2<K,MTag>
  Construct_svd_vertex_2;

  //typedef
  //CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::Construct_sdg_bisector_2<K,MTag>
  //Construct_sdg_bisector_2;

  // Linf traits contain bisector constructions
  typedef Tag_true Tag_has_bisector_constructions;

  template<class Gt, class M>
  class Construct_sdg_bisector_2
   : public CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::
            Construct_sdg_bisector_2<Gt, M>
  {};

  template<class Gt, class M>
  class Construct_sdg_bisector_ray_2
   : public CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::
            Construct_sdg_bisector_ray_2<Gt, M>
  {};

  template<class Gt, class M>
  class Construct_sdg_bisector_segment_2
   : public CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::
            Construct_sdg_bisector_segment_2<Gt, M>
  {};

  // PREDICATES
  //-----------

  // used by triangulation 
  //typedef
  //CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::Orientation_Linf_C2<K>
  //Orientation_2;

  typedef
  CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::Oriented_side_of_bisector_C2<K,MTag>
  Oriented_side_of_bisector_2;

  typedef
  CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::Vertex_conflict_C2<K,MTag>
  Vertex_conflict_2;

  typedef
  CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::Finite_edge_interior_conflict_C2<K,MTag>
  Finite_edge_interior_conflict_2;

  typedef
  CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::Infinite_edge_interior_conflict_C2<K,MTag>
  Infinite_edge_interior_conflict_2;

  typedef
  CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::Oriented_side_C2<K,MTag>
  Oriented_side_2;

public:
  //-----------------------------------------------------------------------
  //                  ACCESS TO OBJECTS
  //-----------------------------------------------------------------------

  // CONSTRUCTIONS
  //--------------
  Construct_svd_vertex_2
  construct_svd_vertex_2_object() const { 
    return Construct_svd_vertex_2();
  }

  //Construct_sdg_bisector_2
  //construct_sdg_bisector_2_object() const {
  //  return Construct_sdg_bisector_2();
  //}

  // PREDICATES
  //-----------

  //Orientation_2
  //orientation_2_object() const {
  //  return Orientation_2();
  //}

  Oriented_side_of_bisector_2
  oriented_side_of_bisector_2_object() const {
    return Oriented_side_of_bisector_2();
  }

  Vertex_conflict_2
  vertex_conflict_2_object() const {
    return Vertex_conflict_2();
  }

  Finite_edge_interior_conflict_2
  finite_edge_interior_conflict_2_object() const {
    return Finite_edge_interior_conflict_2();
  }

  Infinite_edge_interior_conflict_2
  infinite_edge_interior_conflict_2_object() const {
    return Infinite_edge_interior_conflict_2();
  }

  Oriented_side_2
  oriented_side_2_object() const {
    return Oriented_side_2();
  }
};

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_TRAITS_BASE_2_H

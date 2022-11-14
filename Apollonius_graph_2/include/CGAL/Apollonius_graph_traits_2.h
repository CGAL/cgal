// Copyright (c) 2003,2004,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_APOLLONIUS_GRAPH_TRAITS_2_H
#define CGAL_APOLLONIUS_GRAPH_TRAITS_2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Apollonius_graph_2/Predicates_C2.h>
#include <CGAL/Apollonius_graph_2/Constructions_C2.h>
#include <CGAL/Apollonius_graph_2/Incircle8_C2.h>
#ifdef CGAL_APOLLONIUS_GRAPH_D8_TRAITS_2
#include <CGAL/Apollonius_graph_2/Orientation8_C2.h>
#include <CGAL/Apollonius_graph_2/Finite_edge_test8_C2.h>
#endif

#include <CGAL/number_type_basic.h>
#include <CGAL/Apollonius_graph_2/Kernel_wrapper_2.h>

namespace CGAL {

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// the Traits class
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
template < class Rep, class MTag = Integral_domain_without_division_tag >
class Apollonius_graph_traits_2
{
public:
  //-----------------------------------------------------------------------
  //                  TYPE DEFINITIONS
  //-----------------------------------------------------------------------

  // BASIC TYPES
  //------------
private:
  typedef Apollonius_graph_traits_2<Rep,MTag>           Self;
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Apollonius_graph_kernel_wrapper_2<Rep>  Kernel;

public:
  typedef Rep                                           R;
  typedef MTag                                          Method_tag;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Site_2                       Site_2;

  typedef typename Kernel::Line_2                       Line_2;
  typedef typename Kernel::Ray_2                        Ray_2;
  typedef typename Rep::Segment_2                       Segment_2;

  typedef typename Kernel::Object_2                     Object_2;
  typedef typename Kernel::FT                           FT;
  typedef typename Kernel::RT                           RT;

public:
  // OBJECT CONSTRUCTION & ASSIGNMENT
  //---------------------------------
  typedef typename Kernel::Construct_object_2     Construct_object_2;
  typedef typename Kernel::Assign_2               Assign_2;

  // CONSTRUCTIONS
  //--------------
  typedef typename Kernel::Construct_point_2      Construct_point_2;

  // vertex and dual site
  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Construct_Apollonius_vertex_2<Kernel>
  /*                                      */ Construct_Apollonius_vertex_2;

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Construct_Apollonius_site_2<Kernel>
  /*                                        */ Construct_Apollonius_site_2;


  // PREDICATES
  //-----------
  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Compare_x_2<Kernel>   Compare_x_2;

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Compare_y_2<Kernel>   Compare_y_2;

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Compare_weight_2<Kernel>
  Compare_weight_2;

#ifdef CGAL_APOLLONIUS_GRAPH_D8_TRAITS_2
  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Orientation8_C2<Kernel,MTag>
  Orientation_2;
#else
  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Orientation_2<Kernel,MTag>
  Orientation_2;
#endif

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Is_hidden_2<Kernel,MTag>  Is_hidden_2;

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Oriented_side_of_bisector_2<Kernel,MTag>
  /*                                          */ Oriented_side_of_bisector_2;

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Vertex_conflict8_2<Kernel,MTag>
  /*                                                    */ Vertex_conflict_2;

#ifdef CGAL_APOLLONIUS_GRAPH_D8_TRAITS_2
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Finite_edge_interior_conflict8_2<Kernel,MTag>
  /*                                      */ Finite_edge_interior_conflict_2;
#else
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Finite_edge_interior_conflict_2<Kernel,MTag>
  /*                                      */ Finite_edge_interior_conflict_2;
#endif

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Infinite_edge_interior_conflict_2<Kernel,MTag>
  /*                                    */ Infinite_edge_interior_conflict_2;

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Is_degenerate_edge_2<Kernel,MTag>
  /*                                                */  Is_degenerate_edge_2;


public:
  //-----------------------------------------------------------------------
  //                  ACCESS TO OBJECTS
  //-----------------------------------------------------------------------

  // OBJECT CONSTRUCTION & ASSIGNMENT
  Assign_2
  assign_2_object() const {
    return Assign_2();
  }

  Construct_object_2
  construct_object_2_object() const {
    return Construct_object_2();
  }

  Construct_point_2
  construct_point_2_object() const {
    return Construct_point_2();
  }

  // CONSTRUCTIONS
  //--------------
  Construct_Apollonius_vertex_2
  construct_Apollonius_vertex_2_object() const {
    return Construct_Apollonius_vertex_2();
  }

  Construct_Apollonius_site_2
  construct_Apollonius_site_2_object() const {
    return Construct_Apollonius_site_2();
  }


  // PREDICATES
  //-----------
  Compare_x_2
  compare_x_2_object() const {
    return Compare_x_2();
  }

  Compare_y_2
  compare_y_2_object() const {
    return Compare_y_2();
  }

  Compare_weight_2
  compare_weight_2_object() const {
    return Compare_weight_2();
  }

  Orientation_2
  orientation_2_object() const {
    return Orientation_2();
  }

  Is_hidden_2
  is_hidden_2_object() const {
    return Is_hidden_2();
  }

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

  Is_degenerate_edge_2
  is_degenerate_edge_2_object() const {
    return Is_degenerate_edge_2();
  }

};

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_TRAITS_2_H

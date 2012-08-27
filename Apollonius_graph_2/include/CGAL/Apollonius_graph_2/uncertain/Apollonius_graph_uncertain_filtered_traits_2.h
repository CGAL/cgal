// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_APOLLONIUS_GRAPH_UNCERTAIN_FILTERED_TRAITS_2_H
#define CGAL_APOLLONIUS_GRAPH_UNCERTAIN_FILTERED_TRAITS_2_H


#include <CGAL/Apollonius_graph_traits_2.h>
#include <CGAL/Apollonius_graph_2/new_traits/Apollonius_graph_new_traits_2.h>

#include <CGAL/Filtered_predicate.h>
#include <CGAL/Filtered_construction.h>

// includes for the default parameters of the filtered traits
#include <CGAL/Simple_cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/number_utils_classes.h>


// new includes

#include <CGAL/Apollonius_graph_2/uncertain/Uncertain_is_hidden_C2.h>
#include <CGAL/Apollonius_graph_2/uncertain/Uncertain_oriented_side_of_bisector_C2.h>
#include <CGAL/Apollonius_graph_2/uncertain/Uncertain_vertex_conflict_2.h>


namespace CGAL {


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// the filtered Traits class
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

template<class CK_t,
	 class CK_MTag = Integral_domain_without_division_tag,
	 class EK_t    = Simple_cartesian< MP_Float >,
	 class EK_MTag = CK_MTag,
	 class FK_t    = Simple_cartesian< Interval_nt<false> >,
	 class FK_MTag = CK_MTag,
         class C2E_t   = Cartesian_converter<CK_t, EK_t>,
         class C2F_t   =
         Cartesian_converter<CK_t, FK_t, To_interval<typename CK_t::RT> > >
class Apollonius_graph_uncertain_filtered_traits_2
{
private:
  typedef Apollonius_graph_traits_2<CK_t, CK_MTag>    CK_traits;
  typedef Apollonius_graph_traits_2<FK_t, FK_MTag>    FK_traits;
  typedef Apollonius_graph_traits_2<EK_t, EK_MTag>    EK_traits;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Apollonius_graph_kernel_wrapper_2<CK_t>     CK;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Apollonius_graph_kernel_wrapper_2<FK_t>     FK;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Apollonius_graph_kernel_wrapper_2<EK_t>     EK;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::
  Apollonius_graph_cartesian_converter<CK, EK, C2E_t>   C2E;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::
  Apollonius_graph_cartesian_converter<CK, FK, C2F_t>   C2F;

#if 0
  // the following typedefs have been made in the direction of
  // providing filtered constructions; however, there is a problem,
  // namely, the Construct_Apollonius_site_2 functor has two
  // different operator()'s with two different return types; this
  // functor should be split in two (along with the appropriate
  // changes in the spec/concept); see also changes needed for the
  // filtered construction below.
  typedef Cartesian_converter<FK, CK, To_double<typename FK::RT> > F2C_t;
  typedef Cartesian_converter<EK, CK, To_double<typename EK::RT> > E2C_t;

  typedef
  Apollonius_graph_cartesian_converter<FK, CK, F2C_t>   F2C;
  typedef
  Apollonius_graph_cartesian_converter<EK, CK, E2C_t>   E2C;
#endif

  // Types for the construction kernel
  typedef typename CK::Point_2                CK_Point_2;
  typedef typename CK::Site_2                 CK_Site_2;

  typedef typename CK::Line_2                 CK_Line_2;
  typedef typename CK::Ray_2                  CK_Ray_2;
  typedef typename CK::Segment_2              CK_Segment_2;

  typedef typename CK::FT                     CK_FT;
  typedef typename CK::RT                     CK_RT;

  // Types for the exact kernel
  typedef typename EK::Point_2                EK_Point_2;
  typedef typename EK::Site_2                 EK_Site_2;

  typedef typename EK::Line_2                 EK_Line_2;
  typedef typename EK::Ray_2                  EK_Ray_2;
  typedef typename EK::Segment_2              EK_Segment_2;

  typedef typename EK::FT                     EK_FT;
  typedef typename EK::RT                     EK_RT;

  // Types for the filtering kernel
  typedef typename FK::Point_2                FK_Point_2;
  typedef typename FK::Site_2                 FK_Site_2;

  typedef typename FK::Line_2                 FK_Line_2;
  typedef typename FK::Ray_2                  FK_Ray_2;
  typedef typename FK::Segment_2              FK_Segment_2;

  typedef typename FK::FT                     FK_FT;
  typedef typename FK::RT                     FK_RT;

public:
  //-----------------------------------------------------------------------
  //                  TYPE DEFINITIONS
  //-----------------------------------------------------------------------

  // BASIC TYPES
  //------------
  typedef CK_t                          R;
  typedef CK_MTag                       Method_tag;

  typedef CK_traits                     Construction_traits;
  typedef FK_traits                     Filtering_traits;
  typedef EK_traits                     Exact_traits;

  typedef CK_MTag                       Construction_traits_method_tag;
  typedef FK_MTag                       Filtering_traits_method_tag;
  typedef EK_MTag                       Exact_traits_method_tag;

  typedef typename CK::Point_2          Point_2;
  typedef typename CK::Site_2           Site_2;

  typedef typename CK::Line_2           Line_2;
  typedef typename CK::Ray_2            Ray_2;
  typedef typename CK::Segment_2        Segment_2;

  typedef typename CK::Object_2         Object_2;
  typedef typename CK::FT               FT;
  typedef typename CK::RT               RT;


public:
  // OBJECT CONSTRUCTION & ASSIGNMENT
  //---------------------------------
  typedef typename CK_traits::Construct_object_2     Construct_object_2;
  typedef typename CK_traits::Assign_2               Assign_2;

  // CONSTRUCTIONS
  //--------------
  // vertex and dual site
protected:
  typedef typename CK_traits::Construct_Apollonius_vertex_2
  CK_Construct_Apollonius_vertex_2;

  typedef typename CK_traits::Construct_Apollonius_site_2
  CK_Construct_Apollonius_site_2;

  typedef typename FK_traits::Construct_Apollonius_vertex_2
  FK_Construct_Apollonius_vertex_2;

  typedef typename FK_traits::Construct_Apollonius_site_2
  FK_Construct_Apollonius_site_2;

  typedef typename EK_traits::Construct_Apollonius_vertex_2
  EK_Construct_Apollonius_vertex_2;

  typedef typename EK_traits::Construct_Apollonius_site_2
  EK_Construct_Apollonius_site_2;

public:
#if 0
  // the following typedefs have been made in the direction of
  // providing filtered constructions; however, there is a problem,
  // namely, the Construct_Apollonius_site_2 functor has two
  // different operator()'s with two different return types; this
  // functor should be split in two (along with the appropriate
  // changes in the spec/concept); see also changes needed for the
  // filtered construction above.
  typedef Filtered_construction<CK_Construct_Apollonius_vertex_2,
				EK_Construct_Apollonius_vertex_2,
				FK_Construct_Apollonius_vertex_2,
				C2E, C2F, E2C, F2C>
  Construct_Apollonius_vertex_2;

  typedef Filtered_construction<CK_Construct_Apollonius_site_2,
				EK_Construct_Apollonius_site_2,
				FK_Construct_Apollonius_site_2,
				C2E, C2F, E2C, F2C>
  Construct_Apollonius_site_2;
#else
  typedef typename CK_traits::Construct_Apollonius_vertex_2
  Construct_Apollonius_vertex_2;

  typedef typename CK_traits::Construct_Apollonius_site_2
  Construct_Apollonius_site_2;
#endif

private:
  // PREDICATES FOR THE TWO KERNELS
  //-------------------------------

#if 0
  // Predicates for the filtering kernel
  typedef typename FK_traits::Compare_x_2        FK_Compare_x_2;
  typedef typename FK_traits::Compare_y_2        FK_Compare_y_2;
  typedef typename FK_traits::Compare_weight_2   FK_Compare_weight_2;
  typedef typename FK_traits::Orientation_2      FK_Orientation_2;
  typedef typename FK_traits::Is_hidden_2        FK_Is_hidden_2;

  typedef typename FK_traits::Oriented_side_of_bisector_2
  FK_Oriented_side_of_bisector_2;

  typedef typename FK_traits::Vertex_conflict_2  FK_Vertex_conflict_2;

  typedef typename FK_traits::Finite_edge_interior_conflict_2
  FK_Finite_edge_interior_conflict_2;

  typedef typename FK_traits::Infinite_edge_interior_conflict_2
  FK_Infinite_edge_interior_conflict_2;

  typedef typename FK_traits::Is_degenerate_edge_2
  FK_Is_degenerate_edge_2;


  // Predicates for the exact kernel
  typedef typename EK_traits::Compare_x_2        EK_Compare_x_2;
  typedef typename EK_traits::Compare_y_2        EK_Compare_y_2;
  typedef typename EK_traits::Compare_weight_2   EK_Compare_weight_2;
  typedef typename EK_traits::Orientation_2      EK_Orientation_2;
  typedef typename EK_traits::Is_hidden_2        EK_Is_hidden_2;

  typedef typename EK_traits::Oriented_side_of_bisector_2
  EK_Oriented_side_of_bisector_2;

  typedef typename EK_traits::Vertex_conflict_2  EK_Vertex_conflict_2;

  typedef typename EK_traits::Finite_edge_interior_conflict_2
  EK_Finite_edge_interior_conflict_2;

  typedef typename EK_traits::Infinite_edge_interior_conflict_2
  EK_Infinite_edge_interior_conflict_2;

  typedef typename EK_traits::Is_degenerate_edge_2
  EK_Is_degenerate_edge_2;
#else
  // Predicates for the filtering kernel
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Compare_x_2<FK>
  FK_Compare_x_2;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Compare_y_2<FK>
  FK_Compare_y_2;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Compare_weight_2<FK>
  FK_Compare_weight_2;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Orientation_2<FK,FK_MTag>
  FK_Orientation_2;

  //  typedef Ag2_uncertain_is_hidden_C2<FK,FK_MTag> FK_Is_hidden_2;
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Is_hidden_2<FK,FK_MTag>
  FK_Is_hidden_2;

  //  typedef Uncertain_vertex_conflict_new_2<FK,FK_MTag>
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Vertex_conflict_new_2<FK,FK_MTag>
  FK_Vertex_conflict_2;

  //  typedef Ag2_uncertain_oriented_side_of_bisector_C2<FK,FK_MTag>
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Oriented_side_of_bisector_2<FK,FK_MTag>
  FK_Oriented_side_of_bisector_2;

  //  typedef Finite_edge_interior_conflict_new_2<FK,FK_MTag>   
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Finite_edge_interior_conflict_2<FK,FK_MTag>   
  FK_Finite_edge_interior_conflict_2;

  //  typedef Infinite_edge_interior_conflict_new_2<FK,FK_MTag>
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Infinite_edge_interior_conflict_2<FK,FK_MTag>
  FK_Infinite_edge_interior_conflict_2;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Is_degenerate_edge_2<FK,FK_MTag>
  FK_Is_degenerate_edge_2;

  // Predicates for the exact kernel
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Compare_x_2<EK>
  EK_Compare_x_2;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Compare_y_2<EK>
  EK_Compare_y_2;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Compare_weight_2<EK>
  EK_Compare_weight_2;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Orientation_2<EK,EK_MTag>
  EK_Orientation_2;

  //  typedef Ag2_uncertain_is_hidden_C2<EK,EK_MTag> EK_Is_hidden_2;
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Is_hidden_2<EK,EK_MTag>
  EK_Is_hidden_2;

  //  typedef Uncertain_vertex_conflict_new_2<EK,EK_MTag>
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Vertex_conflict_new_2<EK,EK_MTag>
  EK_Vertex_conflict_2;

  //  typedef Ag2_uncertain_oriented_side_of_bisector_C2<EK,EK_MTag>
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Oriented_side_of_bisector_2<EK,EK_MTag>
  EK_Oriented_side_of_bisector_2;

  //  typedef Finite_edge_interior_conflict_new_2<EK,EK_MTag>   
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Finite_edge_interior_conflict_2<EK,EK_MTag>   
  EK_Finite_edge_interior_conflict_2;

  //  typedef Infinite_edge_interior_conflict_new_2<EK,EK_MTag>
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Infinite_edge_interior_conflict_2<EK,EK_MTag>
  EK_Infinite_edge_interior_conflict_2;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Is_degenerate_edge_2<EK,EK_MTag>
  EK_Is_degenerate_edge_2;
#endif

public:
  // PREDICATES
  //-----------

  typedef
  Filtered_predicate<EK_Compare_x_2, FK_Compare_x_2, C2E, C2F>
  Compare_x_2;

  typedef
  Filtered_predicate<EK_Compare_y_2, FK_Compare_y_2, C2E, C2F>
  Compare_y_2;

  typedef
  Filtered_predicate<EK_Compare_weight_2, FK_Compare_weight_2,
		     C2E, C2F>
  Compare_weight_2;

  typedef
  Filtered_predicate<EK_Orientation_2, FK_Orientation_2, C2E, C2F>
  Orientation_2;

  typedef
  Filtered_predicate<EK_Is_hidden_2, FK_Is_hidden_2, C2E, C2F>
  Is_hidden_2;

  typedef
  Filtered_predicate<EK_Oriented_side_of_bisector_2,
			FK_Oriented_side_of_bisector_2, C2E, C2F>
  Oriented_side_of_bisector_2;

  typedef
  Filtered_predicate<EK_Vertex_conflict_2,
			FK_Vertex_conflict_2, C2E, C2F>
  Vertex_conflict_2;

  typedef
  Filtered_predicate<EK_Finite_edge_interior_conflict_2,
		     FK_Finite_edge_interior_conflict_2, C2E, C2F>
  Finite_edge_interior_conflict_2;

  typedef
  Filtered_predicate<EK_Infinite_edge_interior_conflict_2,
		     FK_Infinite_edge_interior_conflict_2, C2E, C2F>
  Infinite_edge_interior_conflict_2;

  typedef
  Filtered_predicate<EK_Is_degenerate_edge_2,
		     FK_Is_degenerate_edge_2, C2E, C2F>
  Is_degenerate_edge_2;

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


#endif // CGAL_APOLLONIUS_GRAPH_UNCERTAIN_FILTERED_TRAITS_2_H

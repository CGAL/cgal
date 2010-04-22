// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_FILTERED_TRAITS_WITH_STATIC_FILTERS_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_FILTERED_TRAITS_WITH_STATIC_FILTERS_2_H

#include <CGAL/Segment_Delaunay_graph_2/basic.h>

#include <CGAL/Segment_Delaunay_graph_2/Filtered_traits_base_2.h>
#include <CGAL/Segment_Delaunay_graph_2/Filtered_traits_concept_check_tags.h>

// includes for the default parameters of the filtered traits
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>
#else
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/number_utils_classes.h>

#include <CGAL/Filtered_predicate.h>
#include <CGAL/Filtered_construction.h>

#include <CGAL/Segment_Delaunay_graph_2/basic.h>
#include <CGAL/Segment_Delaunay_graph_2/Traits_base_2.h>
#ifdef CGAL_SDG_NOX
#include <CGAL/Segment_Delaunay_graph_2/nox/Kernel_wrapper_2.h>
#else
#include <CGAL/Segment_Delaunay_graph_2/Kernel_wrapper_2.h>
#endif
#include <CGAL/Segment_Delaunay_graph_2/Cartesian_converter.h>



CGAL_BEGIN_NAMESPACE

#define SDG2_INS CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Internal
#define SDG2_NS CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// filtered traits class that uses static filters whenever possible
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

// this traits class does NOT support intersecting segments
template<class CK_t  = Simple_cartesian<double>, // the construction kernel
	 class EK_t  = Simple_cartesian<Gmpq>,   // the exact kernel
	 class FK_t  = Simple_cartesian< Interval_nt<false> >,
	 // the filtering kernel
	 class FDK_t  = Filtered_kernel<CK_t>,	 // the filtered kernel
	 class C2E_t = Cartesian_converter<CK_t, EK_t>,
	 class C2F_t =
	 Cartesian_converter<CK_t, FK_t, To_interval<typename CK_t::RT> > >
struct Segment_Delaunay_graph_filtered_traits_with_static_filters_2
{
private:
  typedef
  Segment_Delaunay_graph_filtered_traits_with_static_filters_2
  <CK_t,EK_t,FK_t,FDK_t,C2E_t,C2F_t>
  Self;

  typedef Tag_false                              ITag;
  typedef Field_with_sqrt_tag                    CK_MTag;
  typedef Field_with_sqrt_tag                    FK_MTag;
  typedef Integral_domain_without_division_tag   EK_MTag;

  typedef Segment_Delaunay_graph_traits_base_2<CK_t,CK_MTag,ITag> CK_traits;
  typedef Segment_Delaunay_graph_traits_base_2<FK_t,FK_MTag,ITag> FK_traits;
  typedef Segment_Delaunay_graph_traits_base_2<EK_t,EK_MTag,ITag> EK_traits;

  typedef SDG2_NS::Kernel_wrapper_2<CK_t,ITag>  CK;
  typedef SDG2_NS::Kernel_wrapper_2<FK_t,ITag>  FK;
  typedef SDG2_NS::Kernel_wrapper_2<EK_t,ITag>  EK;

  typedef SDG2_NS::Cartesian_converter<CK,EK,C2E_t>   C2E;
  typedef SDG2_NS::Cartesian_converter<CK,FK,C2F_t>   C2F;

  typedef Cartesian_converter<FK, CK, To_double<typename FK::RT> >  F2C_t;
  typedef SDG2_NS::Cartesian_converter<FK,CK,F2C_t>   F2C;

  typedef Cartesian_converter<EK, CK, To_double<typename EK::RT> >  E2C_t;
  typedef SDG2_NS::Cartesian_converter<EK,CK,E2C_t>   E2C;
  
  // Types for the construction kernel
  typedef typename CK::Point_2                CK_Point_2;
  typedef typename CK::Line_2                 CK_Line_2;
  typedef typename CK::Segment_2              CK_Segment_2;
  typedef typename CK::Ray_2                  CK_Ray_2;

  typedef typename CK::Site_2                 CK_Site_2;

  typedef typename CK::FT                     CK_FT;
  typedef typename CK::RT                     CK_RT;

  // Types for the exact kernel
  typedef typename EK::Point_2                EK_Point_2;
  typedef typename EK::Line_2                 EK_Line_2;
  typedef typename EK::Segment_2              EK_Segment_2;
  typedef typename EK::Ray_2                  EK_Ray_2;

  typedef typename EK::Site_2                 EK_Site_2;

  typedef typename EK::FT                     EK_FT;
  typedef typename EK::RT                     EK_RT;

  // Types for the filtering kernel
  typedef typename FK::Point_2                FK_Point_2;
  typedef typename FK::Line_2                 FK_Line_2;
  typedef typename FK::Segment_2              FK_Segment_2;
  typedef typename FK::Ray_2                  FK_Ray_2;

  typedef typename FK::Site_2                 FK_Site_2;

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
  typedef ITag                          Intersections_tag;

  typedef CK_t                          Construction_kernel;
  typedef FK_t                          Filtering_kernel;
  typedef FDK_t                         Filtered_kernel;
  typedef EK_t                          Exact_kernel;

  typedef FDK_t                         FDK;

  typedef CK_traits                     Construction_traits;
  typedef FK_traits                     Filtering_traits;
  typedef EK_traits                     Exact_traits;

  //  typedef CK_MTag                       Construction_traits_method_tag;
  //  typedef FK_MTag                       Filtering_traits_method_tag;
  //  typedef EK_MTag                       Exact_traits_method_tag;

  typedef typename CK::Point_2          Point_2;
  typedef typename CK::Line_2           Line_2;
  typedef typename CK::Segment_2        Segment_2;
  typedef typename CK::Ray_2            Ray_2;
  //  typedef typename CK::Circle_2         Circle_2;
  typedef typename CK::Site_2           Site_2;

  typedef typename CK::Object_2         Object_2;

  typedef typename CK::FT               FT;
  typedef typename CK::RT               RT;

  typedef typename CK::Rep_tag          Rep_tag;

protected:
  typedef
  CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Internal::Arrangement_enum
  Arrangement_enum;

private:
  typedef typename CK_traits::Construct_svd_vertex_2
  CK_Construct_svd_vertex_2;

  typedef typename FK_traits::Construct_svd_vertex_2
  FK_Construct_svd_vertex_2;

  typedef typename EK_traits::Construct_svd_vertex_2
  EK_Construct_svd_vertex_2;

public:
  // OBJECT CONSTRUCTION & ASSIGNMENT
  //---------------------------------
  typedef typename CK::Construct_object_2     Construct_object_2;
  typedef typename CK::Assign_2               Assign_2;

  // CONSTRUCTIONS
  //--------------
  // vertex and Voronoi circle
  typedef
  Filtered_construction<CK_Construct_svd_vertex_2,
			EK_Construct_svd_vertex_2,
			FK_Construct_svd_vertex_2,
			C2E, C2F, E2C, F2C>
  Construct_svd_vertex_2;

private:
  // PREDICATES FOR THE TWO KERNELS
  //-------------------------------

  // Predicates for the filtering kernel
  typedef typename FK_traits::Compare_x_2        FK_Compare_x_2;
  typedef typename FK_traits::Compare_y_2        FK_Compare_y_2;
  typedef typename FK_traits::Orientation_2      FK_Orientation_2;
  typedef typename FK_traits::Equal_2            FK_Equal_2;
  typedef typename FK_traits::Are_parallel_2     FK_Are_parallel_2;

  typedef typename FK_traits::Oriented_side_of_bisector_2
  FK_Oriented_side_of_bisector_2;

  typedef typename FK_traits::Vertex_conflict_2  FK_Vertex_conflict_2;

  typedef typename FK_traits::Finite_edge_interior_conflict_2
  FK_Finite_edge_interior_conflict_2;

  typedef typename FK_traits::Infinite_edge_interior_conflict_2
  FK_Infinite_edge_interior_conflict_2;

  typedef typename FK_traits::Is_degenerate_edge_2
  FK_Is_degenerate_edge_2;

  typedef typename FK_traits::Arrangement_type_2 FK_Arrangement_type_2;
  typedef typename FK_traits::Oriented_side_2    FK_Oriented_side_2;

  // Predicates for the exact kernel
  typedef typename EK_traits::Compare_x_2        EK_Compare_x_2;
  typedef typename EK_traits::Compare_y_2        EK_Compare_y_2;
  typedef typename EK_traits::Orientation_2      EK_Orientation_2;
  typedef typename EK_traits::Equal_2            EK_Equal_2;
  typedef typename EK_traits::Are_parallel_2     EK_Are_parallel_2;

  typedef typename EK_traits::Oriented_side_of_bisector_2
  EK_Oriented_side_of_bisector_2;

  typedef typename EK_traits::Vertex_conflict_2  EK_Vertex_conflict_2;

  typedef typename EK_traits::Finite_edge_interior_conflict_2
  EK_Finite_edge_interior_conflict_2;

  typedef typename EK_traits::Infinite_edge_interior_conflict_2
  EK_Infinite_edge_interior_conflict_2;

  typedef typename EK_traits::Is_degenerate_edge_2
  EK_Is_degenerate_edge_2;

  typedef typename EK_traits::Arrangement_type_2 EK_Arrangement_type_2;
  typedef typename EK_traits::Oriented_side_2    EK_Oriented_side_2;


public:
  // PREDICATES
  //-----------

  struct C2FD
  {
    inline
    typename FDK::Point_2 operator()(const Point_2& p) const
    {
      return typename FDK::Point_2(p.x(), p.y());
    }
  };

#if 1
  struct Compare_x_2
  {
    inline Comparison_result
    operator()(const Site_2& s, const Site_2& t) const
    {
      CGAL_precondition( s.is_point() && t.is_point() );
      C2FD c2fd;
      return typename FDK::Compare_x_2()(c2fd(s.point()), c2fd(t.point()));
    }

    inline Comparison_result
    operator()(const Point_2& p, const Point_2& q) const
    {
      C2FD c2fd;
      return typename FDK::Compare_x_2()(c2fd(p), c2fd(q));
    }
  };
#else
  typedef
  Filtered_predicate<EK_Compare_x_2, FK_Compare_x_2, C2E, C2F>
  Compare_x_2;
#endif

#if 1
  struct Compare_y_2
  {
    inline Comparison_result 
    operator()(const Site_2& s, const Site_2& t) const
    {
      CGAL_precondition( s.is_point() && t.is_point() );
      C2FD c2fd;
      return typename FDK::Compare_y_2()(c2fd(s.point()), c2fd(t.point()));
    }

    inline Comparison_result
    operator()(const Point_2& p, const Point_2& q) const
    {
      C2FD c2fd;
      return typename FDK::Compare_y_2()(c2fd(p), c2fd(q));
    }
  };
#else
  typedef
  Filtered_predicate<EK_Compare_y_2, FK_Compare_y_2, C2E, C2F>
  Compare_y_2;
#endif

#if 1
  struct Orientation_2
  {
    inline Orientation
    operator()(const Site_2& s1, const Site_2& s2, const Site_2& s3) const
    {
      CGAL_precondition( s1.is_point() && s2.is_point() && s3.is_point() );
      C2FD c2fd;
      return typename FDK::Orientation_2()(c2fd(s1.point()),
					   c2fd(s2.point()),
					   c2fd(s3.point()));
    }

    inline Orientation
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      C2FD c2fd;
      return typename FDK::Orientation_2()(c2fd(p),
					   c2fd(q),
					   c2fd(r));
    }
  };
#else
  typedef
  Filtered_predicate<EK_Orientation_2, FK_Orientation_2, C2E, C2F>
  Orientation_2;
#endif

#if 1
  struct Equal_2
  {
    inline bool
    operator()(const Site_2& s, const Site_2& t) const
    {
      CGAL_precondition( s.is_point() && t.is_point() );
      C2FD c2fd;
      return typename FDK::Equal_2()(c2fd(s.point()), c2fd(t.point()));
    }

    inline bool
    operator()(const Point_2& p, const Point_2& q) const
    {
      C2FD c2fd;
      return typename FDK::Equal_2()(c2fd(p), c2fd(q));
    }
  };
#else
  typedef
  Filtered_predicate<EK_Equal_2, FK_Equal_2, C2E, C2F>
  Equal_2;
#endif

#if 0
  struct Are_parallel_2
  {
    inline bool
    operator()(const Site_2& s, const Site_2& t) const
    {
      CGAL_precondition( s.is_segment() && t.is_segment() );
      return (FDK::Compare_slope_2()(s.segment(), t.segment()) == EQUAL);
    }
  };
#else
  typedef
  Filtered_predicate<EK_Are_parallel_2,	FK_Are_parallel_2, C2E, C2F>
  Are_parallel_2;
#endif

#if 1
  class Oriented_side_of_bisector_2
  {
  private:
    typedef typename FK_traits::Oriented_side_of_bisector_2 Filtering_predicate;
    typedef typename EK_traits::Oriented_side_of_bisector_2 Exact_predicate;

    typedef Filtered_predicate<Exact_predicate,Filtering_predicate,C2E,C2F>
    Filtered_exact_predicate;

  public:
    inline Oriented_side
    operator()(const Site_2& s1, const Site_2& s2, const Site_2& p) const
    {
      C2FD c2fd;
      if ( p.is_point() && s1.is_point() && s2.is_point() ) {
	Comparison_result cr =
	  typename FDK::Compare_distance_2()(c2fd(p.point()),
					     c2fd(s1.point()),
					     c2fd(s2.point()));
	if ( cr == EQUAL ) { return ON_ORIENTED_BOUNDARY; }
	if ( cr == SMALLER ) { return ON_POSITIVE_SIDE; }
	return ON_NEGATIVE_SIDE;
      }
      return Filtered_exact_predicate()(s1, s2, p);
    }
  };
#else
  typedef
  Filtered_predicate<EK_Oriented_side_of_bisector_2,
  		     FK_Oriented_side_of_bisector_2, C2E, C2F>
  Oriented_side_of_bisector_2;
#endif

#if 1
  class Vertex_conflict_2
  {
  private:
    typedef typename FK_traits::Vertex_conflict_2   Filtering_predicate;
    typedef typename EK_traits::Vertex_conflict_2   Exact_predicate;

    typedef Filtered_predicate<Exact_predicate,Filtering_predicate,C2E,C2F>
    Filtered_exact_predicate;

  public:
    inline Sign
    operator()(const Site_2& s1, const Site_2& s2,
	       const Site_2& s3, const Site_2& q) const
    {
      if ( s1.is_point() && s2.is_point() && s3.is_point() && q.is_point() ) {
	C2FD c2fd;
	return -(typename FDK::Side_of_oriented_circle_2()(c2fd(s1.point()),
							   c2fd(s2.point()),
							   c2fd(s3.point()),
							   c2fd(q.point())));
      }
      return Filtered_exact_predicate()(s1, s2, s3, q);
    }

    inline Sign
    operator()(const Point_2& p1, const Point_2& p2,
	       const Point_2& p3, const Point_2& q) const
    {
      C2FD c2fd;
      return -(typename FDK::Side_of_oriented_circle_2()(c2fd(p1),
							 c2fd(p2),
							 c2fd(p3),
							 c2fd(q)));
  }

    inline Sign
    operator()(const Site_2& s1, const Site_2& s2,
	       const Site_2& q) const
    {
      return Filtered_exact_predicate()(s1, s2, q);
    }

  };
#else
  typedef
  Filtered_predicate<EK_Vertex_conflict_2,
  		     FK_Vertex_conflict_2, C2E, C2F>
  Vertex_conflict_2;
#endif

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

private:
  typedef
  Filtered_predicate<EK_Arrangement_type_2, FK_Arrangement_type_2, C2E, C2F>
  Arrangement_type_2_base;

public:
  struct Arrangement_type_2
    : public Arrangement_type_2_base, public Arrangement_enum
  {};

  typedef
  Filtered_predicate<EK_Oriented_side_2, FK_Oriented_side_2, C2E, C2F>
  Oriented_side_2;

public:
  //-----------------------------------------------------------------------
  //                  ACCESS TO OBJECTS
  //-----------------------------------------------------------------------

  // OBJECT CONSTRUCTION & ASSIGNMENT
  //---------------------------------
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
  Construct_svd_vertex_2
  construct_svd_vertex_2_object() const { 
    return Construct_svd_vertex_2();
  }

  /*
  Construct_site_2
  construct_site_2_object() const { 
    return Construct_site_2();
  }
  */

  /*
  Construct_sdg_circle_2
  construct_sdg_circle_2_object() const {
    return Construct_sdg_circle_2();
  }
  */

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

  Orientation_2
  orientation_2_object() const {
    return Orientation_2();
  }

  Equal_2
  equal_2_object() const {
    return Equal_2();
  }

  Are_parallel_2
  are_parallel_2_object() const {
    return Are_parallel_2();
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

  Arrangement_type_2
  arrangement_type_2_object() const {
    return Arrangement_type_2();
  }

  Oriented_side_2
  oriented_side_2_object() const {
    return Oriented_side_2();
  }
};

#undef SDG2_NS

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_FILTERED_TRAITS_WITH_STATIC_FILTERS_2_H

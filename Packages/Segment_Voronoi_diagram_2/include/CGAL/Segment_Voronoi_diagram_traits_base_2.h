// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_BASE_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_BASE_2_H

#include <CGAL/Segment_Voronoi_diagram_short_names_2.h>

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/predicates/Segment_Voronoi_diagram_predicates_C2.h>
#include <CGAL/predicates/Segment_Voronoi_diagram_predicates_ftC2.h>
#include <CGAL/Segment_Voronoi_diagram_constructions_C2.h>
#endif


#include <CGAL/Number_type_traits.h>
#include <CGAL/Segment_Voronoi_diagram_kernel_wrapper_2.h>




CGAL_BEGIN_NAMESPACE

//***********************************************************************
//***********************************************************************
//                              PREDICATES
//***********************************************************************
//***********************************************************************

//-----------------------------------------------------------------------
//                           compare x
//-----------------------------------------------------------------------

template< class K >
class Svd_compare_x_2
{
public:
  typedef typename K::Site_2         Site_2;
  typedef typename K::Point_2        Point_2;

  typedef Comparison_result          result_type;

  typedef Arity_tag<2>               Arity;

private:
  typedef typename K::Compare_x_2    compare_x_2;

public:

  result_type operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );
    return compare_x_2()( p.point(), q.point() );
  }
};

//-----------------------------------------------------------------------
//                           compare y
//-----------------------------------------------------------------------

template< class K >
class Svd_compare_y_2
{
public:
  typedef typename K::Site_2         Site_2;
  typedef typename K::Point_2        Point_2;
  typedef Comparison_result          result_type;
  typedef Arity_tag<2>               Arity;

private:
  typedef typename K::Compare_y_2  compare_y_2;

public:

  result_type operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );
    return compare_y_2()( p.point(), q.point() );
  }
};

//-----------------------------------------------------------------------
//                           orientation
//-----------------------------------------------------------------------

template< class K >
class Svd_orientation_2
{
public:
  typedef typename K::Site_2         Site_2;
  typedef typename K::Point_2        Point_2;
  typedef Orientation                result_type;
  typedef Arity_tag<3>               Arity;

private:
  typedef typename K::Orientation_2  orientation_2;

public:

  result_type operator()(const Site_2& p, const Site_2& q,
			 const Site_2& r) const
  {
    CGAL_precondition( p.is_point() && q.is_point() && r.is_point() );
    return svd_orientation_C2<K>( p, q, r );
  }
};


//-----------------------------------------------------------------------
//                           are same points
//-----------------------------------------------------------------------

template< class K >
class Are_same_points_2
{
public:
  typedef typename K::Site_2        Site_2;
  typedef typename K::Point_2       Point_2;
  typedef bool                      result_type;
  typedef Arity_tag<2>              Arity;

public:

  result_type operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );
    return svd_are_same_points_C2<K>(p, q);
  }
};

//-----------------------------------------------------------------------
//                           are parallel
//-----------------------------------------------------------------------


template< class K >
class Are_parallel_2
{
public:
  typedef typename K::Site_2       Site_2;
  typedef bool                     result_type;
  typedef Arity_tag<2>             Arity;

public:
  result_type operator()(const Site_2& p, const Site_2& q) const
  {
    return svd_are_parallel_C2(p, q);
  }
};



//-----------------------------------------------------------------------
//                       oriented side of bisector
//-----------------------------------------------------------------------

template<class K, class Method_tag>
class Svd_oriented_side_of_bisector_2
{
public:
  typedef typename K::Site_2     Site_2;
  typedef Oriented_side          result_type;
  typedef Arity_tag<3>           Arity;

public:
  result_type operator()(const Site_2& t1, const Site_2& t2,
			 const Site_2& q) const
  {
    return svd_oriented_side_of_bisector_ftC2<K,Method_tag>
      (t1, t2, q, Method_tag());
  }

};

//-----------------------------------------------------------------------
//                           vertex conflict
//-----------------------------------------------------------------------


template< class K, class Method_tag >
class Svd_vertex_conflict_2
{
public:
  typedef typename K::Site_2  Site_2;
  typedef Sign                result_type;

  struct Arity {};

public:
  result_type operator()(const Site_2& p, const Site_2& q,
			 const Site_2& t) const
  {
    return
      svd_vertex_conflict_ftC2<K,Method_tag>(p, q, t, Method_tag());
  }

  result_type operator()(const Site_2& p, const Site_2& q,
			 const Site_2& r, const Site_2& t) const
  {
    return
      svd_vertex_conflict_ftC2<K,Method_tag>(p, q, r, t, Method_tag());
  }
};

//-----------------------------------------------------------------------
//                   finite edge interior conflict
//-----------------------------------------------------------------------

template< class K, class Method_tag >
class Svd_finite_edge_interior_conflict_2
{
public:
  typedef typename K::Site_2  Site_2;
  typedef bool                result_type;

  struct Arity {};

public:
  result_type operator()(const Site_2& p, const Site_2& q,
			 const Site_2& t, Sign sgn) const
  {
    return
      svd_finite_edge_conflict_ftC2<K,Method_tag>(p, q, t,
						  sgn, Method_tag());
  }


  result_type
  operator()(const Site_2& p, const Site_2& q,
	     const Site_2& r, const Site_2& t, Sign sgn) const
  {
    return
      svd_finite_edge_conflict_ftC2<K,Method_tag>(p, q, r, t,
						  sgn, Method_tag());
  }

  result_type operator()(const Site_2& p, const Site_2& q,
			 const Site_2& r, const Site_2& s,
			 const Site_2& t, Sign sgn) const
  {
    return
      svd_finite_edge_conflict_ftC2<K,Method_tag>(p, q, r, s, t,
						  sgn, Method_tag());
  }

};


//-----------------------------------------------------------------------
//                   infinite edge interior conflict
//-----------------------------------------------------------------------


template< class K, class Method_tag >
class Svd_infinite_edge_interior_conflict_2
{
public:
  typedef typename K::Site_2  Site_2;
  typedef bool                result_type;
  typedef Arity_tag<5>        Arity;

public:
  result_type
  operator()(const Site_2& q, const Site_2& r, const Site_2& s,
	     const Site_2& t, Sign sgn) const
  {
    return
      svd_infinite_edge_conflict_ftC2<K,Method_tag>(q, r, s, t,
						    sgn, Method_tag());
  }

};

//-----------------------------------------------------------------------
//                          is degenerate edge
//-----------------------------------------------------------------------

template< class K, class Method_tag >
class Svd_is_degenerate_edge_2
{
public:
  typedef typename K::Site_2  Site_2;
  typedef bool                result_type;
  typedef Arity_tag<4>        Arity;

public:
  result_type operator()(const Site_2& p, const Site_2& q,
			 const Site_2& r, const Site_2& s) const
  {
    return
      svd_is_degenerate_edge_ftC2<K,Method_tag>(p, q, r, s, Method_tag());
  }
};

//-----------------------------------------------------------------------
//                          do intersect
//-----------------------------------------------------------------------

template< class K, class Method_tag >
class Svd_arrangement_type_2
{
public:
  typedef typename K::Site_2     Site_2;
  typedef std::pair<int,int>     result_type;
  typedef Arity_tag<2>           Arity;

public:
  result_type operator()(const Site_2& p, const Site_2& q) const
  {
    return
      svd_arrangement_type_C2<K,Method_tag>(p, q, Method_tag());
  }
};


//-----------------------------------------------------------------------
//                          oriented side
//-----------------------------------------------------------------------

template< class K, class Method_tag >
class Svd_oriented_side_2
{
public:
  typedef typename K::Site_2     Site_2;
  typedef Oriented_side          result_type;
  typedef Arity_tag<5>           Arity;

public:
  result_type
  operator()(const Site_2& s1, const Site_2& s2, const Site_2& s3,
	     const Site_2& s, const Site_2& p) const
  {
    return
      svd_oriented_side_ftC2<K,Method_tag>(s1, s2, s3, s, p, Method_tag());
  }
};



//-----------------------------------------------------------------------
// the Traits class
//-----------------------------------------------------------------------
template<class R, class MTag, class ITag>
class Segment_Voronoi_diagram_traits_base_2
{
public:
  //-----------------------------------------------------------------------
  //                  TYPE DEFINITIONS
  //-----------------------------------------------------------------------

  // BASIC TYPES
  //------------
  
  typedef Segment_Voronoi_diagram_kernel_wrapper_2<R,ITag> Kernel;
  typedef Kernel                                           K;
  typedef R                                                Rep;
  typedef MTag                                             Method_tag;

  // the following tag controls how support intersections in the
  // traits. If it is Tag_true then we fully support intersections.
  // If it is Tag_true it is assumed that no intersections appear in
  // the data and so there is limited support for intersections.
  typedef ITag                                    Intersections_tag;

  typedef typename Kernel::Point_2                Point_2;
  typedef typename Kernel::Line_2                 Line_2;
  typedef typename Kernel::Segment_2              Segment_2;
  typedef typename Kernel::Ray_2                  Ray_2;
  //  typedef typename Kernel::Circle_2               Circle_2;

  typedef typename Kernel::Site_2                 Site_2;
  typedef typename Kernel::Object_2               Object_2;

  typedef typename Kernel::FT                     FT;
  typedef typename Kernel::RT                     RT;

public:
// OBJECT CONSTRUCTION & ASSIGNMENT
  //-------------------------------
  typedef typename Kernel::Construct_object_2     Construct_object_2;
  typedef typename Kernel::Assign_2               Assign_2;

  // CONSTRUCTIONS
  //--------------
  // vertex and Voronoi circle
  typedef CGAL::Construct_svd_vertex_2<K,MTag>  Construct_svd_vertex_2;
  //  typedef CGAL::Construct_svd_circle_2<K,MTag>  Construct_svd_circle_2;

  // PREDICATES
  //-----------
  typedef CGAL::Svd_compare_x_2<K>                      Compare_x_2;
  typedef CGAL::Svd_compare_y_2<K>                      Compare_y_2;
  typedef CGAL::Svd_orientation_2<K>                    Orientation_2;
  typedef CGAL::Are_same_points_2<K>                    Are_same_points_2;
  typedef CGAL::Are_parallel_2<K>                       Are_parallel_2;
  typedef CGAL::Svd_oriented_side_of_bisector_2<K,MTag> 
  /*                                          */ Oriented_side_of_bisector_2;
  typedef CGAL::Svd_vertex_conflict_2<K,MTag >             Vertex_conflict_2;
  typedef CGAL::Svd_finite_edge_interior_conflict_2<K,MTag >
  /*                                      */ Finite_edge_interior_conflict_2;
  typedef CGAL::Svd_infinite_edge_interior_conflict_2<K,MTag>
  /*                                    */ Infinite_edge_interior_conflict_2;
  typedef CGAL::Svd_is_degenerate_edge_2<K,MTag>        Is_degenerate_edge_2;
  typedef CGAL::Svd_arrangement_type_2<K,MTag>          Arrangement_type_2;
  typedef CGAL::Svd_oriented_side_2<K,MTag>             Oriented_side_2;


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
  Construct_svd_circle_2
  construct_svd_circle_2_object() const {
    return Construct_svd_circle_2();
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

  Are_same_points_2
  are_same_points_2_object() const {
    return Are_same_points_2();
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


CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_BASE_2_H

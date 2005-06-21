// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
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



#ifndef CGAL_APOLLONIUS_GRAPH_TRAITS_2_H
#define CGAL_APOLLONIUS_GRAPH_TRAITS_2_H


#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/predicates/Apollonius_graph_ftC2.h>
#endif

#include <CGAL/Number_type_traits.h>
#include <CGAL/Apollonius_graph_kernel_wrapper_2.h>



CGAL_BEGIN_NAMESPACE


//***********************************************************************
//***********************************************************************
//                              PREDICATES
//***********************************************************************
//***********************************************************************

//-----------------------------------------------------------------------
//                        Orientation
//-----------------------------------------------------------------------

template< class K >
class AG2_Orientation_test_2
{
 public:
  typedef typename K::Site_2   Site_2;
  typedef Orientation          result_type;
  typedef Arity_tag<3>         Arity;

  inline Orientation operator()(const Site_2 &p, const Site_2 &q,
				const Site_2 &r) const
  {
    return ag2_orientation_test_C2(p.x(), p.y(), p.weight(),
				   q.x(), q.y(), q.weight(),
				   r.x(), r.y(), r.weight());
  }
};

//-----------------------------------------------------------------------
//                        Is hidden
//-----------------------------------------------------------------------

template< class K, class Method_tag >
class Is_hidden_2
{
 public:
  typedef typename K::Site_2   Site_2;
  typedef bool                 result_type;
  typedef Arity_tag<2>         Arity;

 private:
  bool is_hidden(const Site_2 &p, const Site_2 &q, Ring_tag) const
  {
    return ag_is_hidden_test_ring_C2(p.x(), p.y(), p.weight(),
				     q.x(), q.y(), q.weight());

  }

  bool is_hidden(const Site_2 &p, const Site_2 &q, Sqrt_field_tag) const
  {
    return ag_is_hidden_test_sqrtf_C2(p.x(), p.y(), p.weight(),
				      q.x(), q.y(), q.weight());
  }

 public:
  inline bool operator()(const Site_2 &p, const Site_2 &q) const
  {
    return is_hidden(p, q, Method_tag());
  }
};


//-----------------------------------------------------------------------
//                    Oriented side of bisector
//-----------------------------------------------------------------------


template< class K, class Method_tag >
class Oriented_side_of_bisector_2
{
 public:
  typedef typename K::Point_2             Point_2;
  typedef typename K::Site_2              Site_2;
  typedef Oriented_side                   result_type;
  typedef Arity_tag<3>                    Arity;

 private:
  Comparison_result compare_distances(const Site_2& p1,
				      const Site_2& p2,
				      const Point_2 &p,
				      Ring_tag) const
  {
    return compare_ag_distances_test_ring_C2(p1.x(), p1.y(), p1.weight(),
					     p2.x(), p2.y(), p2.weight(),
					     p.x(),  p.y());
  }

  Comparison_result compare_distances(const Site_2& p1,
				      const Site_2& p2,
				      const Point_2 &p,
				      Sqrt_field_tag) const
  {
    return
      compare_ag_distances_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
					 p2.x(), p2.y(), p2.weight(),
					 p.x(),  p.y());
  }

 public:
  inline Oriented_side operator()(const Site_2& p1,
				  const Site_2& p2,
				  const Point_2 &p) const
  {
    Comparison_result r = compare_distances(p1, p2, p, Method_tag());

    if ( r == EQUAL ) { return ON_ORIENTED_BOUNDARY; }
    return ( r == LARGER ) ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE;
  }
};



//-----------------------------------------------------------------------
//                        Vertex conflict
//-----------------------------------------------------------------------

template < class K, class Method_tag >
class Vertex_conflict_2
{
 public:
  typedef typename K::Site_2      Site_2;
  typedef Sign                    result_type;
  struct Arity {};

 private:
  Sign incircle(const Site_2& p1, const Site_2& p2,
		const Site_2& p3, const Site_2& q, 
		Sqrt_field_tag) const
  {
    return ag_incircle_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
				     p2.x(), p2.y(), p2.weight(),
				     p3.x(), p3.y(), p3.weight(),
				      q.x(),  q.y(),  q.weight());
  }

  Sign incircle(const Site_2& p1, const Site_2& p2,
		const Site_2& p3, const Site_2& q, 
		Ring_tag) const
  {
    return ag_incircle_test_ring_C2(p1.x(), p1.y(), p1.weight(),
				    p2.x(), p2.y(), p2.weight(),
				    p3.x(), p3.y(), p3.weight(),
				     q.x(),  q.y(),  q.weight());
  }

  Sign incircle(const Site_2& p1, const Site_2& p2,
		const Site_2& q,
		Sqrt_field_tag) const
  {
    return ag_incircle_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
				     p2.x(), p2.y(), p2.weight(),
				      q.x(),  q.y(),  q.weight());
  }

  Sign incircle(const Site_2& p1, const Site_2& p2,
		const Site_2& q, 
		Ring_tag) const
  {
    return ag_incircle_test_ring_C2(p1.x(), p1.y(), p1.weight(),
				    p2.x(), p2.y(), p2.weight(),
				     q.x(),  q.y(),  q.weight());
  }

 public:
  inline
  Sign operator()(const Site_2& p1, const Site_2& p2,
		  const Site_2& p3, const Site_2& q) const
  {
    typedef typename K::Rep_tag Tag;
    return incircle(p1, p2, p3, q, Method_tag());
  }


  inline
  Sign operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& q) const
  {
    typedef typename K::Rep_tag Tag;
    return incircle(p1, p2, q, Method_tag());
  }
};

//-----------------------------------------------------------------------
//                    Finite edge interior conflict
//-----------------------------------------------------------------------

template < class K, class Method_tag >
class Finite_edge_interior_conflict_2
{
 public:
  typedef typename K::Site_2  Site_2;
  typedef bool                result_type;
  struct Arity {};

 private:
  bool finite_edge_test(const Site_2& p1, const Site_2& p2,
			const Site_2& p3, const Site_2& p4,
			const Site_2& q, bool b,
			Sqrt_field_tag) const
  {
  return ag_finite_edge_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
				      p2.x(), p2.y(), p2.weight(),
				      p3.x(), p3.y(), p3.weight(),
				      p4.x(), p4.y(), p4.weight(),
				       q.x(),  q.y(), q.weight(), b);
  }

  bool finite_edge_test(const Site_2& p1, const Site_2& p2,
			const Site_2& p3, const Site_2& p4,
			const Site_2& q, bool b,
			Ring_tag) const
  {
    return ag_finite_edge_test_ring_C2(p1.x(), p1.y(), p1.weight(),
				       p2.x(), p2.y(), p2.weight(),
				       p3.x(), p3.y(), p3.weight(),
				       p4.x(), p4.y(), p4.weight(),
				        q.x(),  q.y(), q.weight(), b);
  }

  bool finite_edge_test(const Site_2& p1, const Site_2& p2,
			const Site_2& p3, const Site_2& q, bool b,
			Sqrt_field_tag) const
  {
    return ag_finite_edge_test_degenerated_sqrtf_C2(p1.x(), p1.y(),
						    p1.weight(),
						    p2.x(), p2.y(),
						    p2.weight(),
						    p3.x(), p3.y(),
						    p3.weight(),
						    q.x(),  q.y(),
						    q.weight(), b);
  }

  bool finite_edge_test(const Site_2& p1, const Site_2& p2,
			const Site_2& p3, const Site_2& q, bool b,
			Ring_tag) const
  {
    return ag_finite_edge_test_degenerated_ring_C2(p1.x(), p1.y(),
						   p1.weight(),
						   p2.x(), p2.y(),
						   p2.weight(),
						   p3.x(), p3.y(),
						   p3.weight(),
						   q.x(),  q.y(),
						   q.weight(), b);
  }

  bool finite_edge_test(const Site_2& p1, const Site_2& p2,
			const Site_2& q, bool b,
			Sqrt_field_tag) const
  {
    return ag_finite_edge_test_degenerated_sqrtf_C2(p1.x(), p1.y(),
						    p1.weight(),
						    p2.x(), p2.y(),
						    p2.weight(),
						    q.x(),  q.y(),
						    q.weight(), b);
  }

  bool finite_edge_test(const Site_2& p1, const Site_2& p2,
			const Site_2& q, bool b,
			Ring_tag) const
  {
    return ag_finite_edge_test_degenerated_ring_C2(p1.x(), p1.y(),
						   p1.weight(),
						   p2.x(), p2.y(),
						   p2.weight(),
						   q.x(),  q.y(),
						   q.weight(), b);
  }
  
 public:
  inline
  bool operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& q, bool b) const
  {
    return finite_edge_test(p1, p2, p3, q, b, Method_tag());
  }

  bool operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& q, bool b) const
  {
    return finite_edge_test(p1, p2, q, b, Method_tag());
  }



  inline
  bool operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& p4,
		  const Site_2& q,
		  bool b) const
  {
    return finite_edge_test(p1, p2, p3, p4, q, b, Method_tag());
  }
};

//-----------------------------------------------------------------------
//                   Infinite edge interior conflict
//-----------------------------------------------------------------------


template < class K, class Method_tag >
class Infinite_edge_interior_conflict_2
{
 public:
  typedef typename K::Site_2   Site_2;
  typedef bool                 result_type;
  typedef Arity_tag<5>         Arity;

 private:
  bool infinite_edge_test(const Site_2& p2, const Site_2& p3,
			  const Site_2& p4, const Site_2& q, bool b,
			  Sqrt_field_tag) const
  {
    return
      ag_infinite_edge_test_sqrtf_C2(p2.x(), p2.y(), p2.weight(),
				     p3.x(), p3.y(), p3.weight(),
				     p4.x(), p4.y(), p4.weight(),
				      q.x(),  q.y(),  q.weight(), b);
  }

  bool infinite_edge_test(const Site_2& p2, const Site_2& p3,
			  const Site_2& p4, const Site_2& q, bool b,
			  Ring_tag) const
  {
    return
      ag_infinite_edge_test_ring_C2(p2.x(), p2.y(), p2.weight(),
				    p3.x(), p3.y(), p3.weight(),
				    p4.x(), p4.y(), p4.weight(),
				     q.x(),  q.y(),  q.weight(), b);
  }

 public:
  bool operator()(const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& p4,
		  const Site_2& q, bool b) const
  {
    return infinite_edge_test(p2, p3, p4, q, b, Method_tag());
  }
};


//-----------------------------------------------------------------------
//                          Is degenerate
//-----------------------------------------------------------------------


template < class K, class Method_tag >
class Is_degenerate_edge_2
{
 public:
  typedef typename K::Site_2    Site_2;
  typedef bool                  result_type;
  typedef Arity_tag<4>          Arity;

 private:
  bool is_degenerate_edge(const Site_2& p1, const Site_2& p2,
			  const Site_2& p3, const Site_2& p4,
			  Sqrt_field_tag) const
  {
    return
      ag_is_degenerate_edge_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
					  p2.x(), p2.y(), p2.weight(),
					  p3.x(), p3.y(), p3.weight(),
					  p4.x(), p4.y(), p4.weight());
  }

  bool is_degenerate_edge(const Site_2& p1, const Site_2& p2,
			  const Site_2& p3, const Site_2& p4,
			  Ring_tag) const
  {
    return
      ag_is_degenerate_edge_test_ring_C2(p1.x(), p1.y(), p1.weight(),
					 p2.x(), p2.y(), p2.weight(),
					 p3.x(), p3.y(), p3.weight(),
					 p4.x(), p4.y(), p4.weight());
  }

 public:
  inline bool operator()(const Site_2& p1,
			 const Site_2& p2,
			 const Site_2& p3,
			 const Site_2& p4) const
  {
    return is_degenerate_edge(p1, p2, p3, p4, Method_tag());
  }
};


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// the Traits class
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
template < class Rep, class MTag = Ring_tag >
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
  typedef Apollonius_graph_kernel_wrapper_2<Rep>        Kernel;

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
  // vertex and dual site
  typedef CGAL::Construct_Apollonius_vertex_2<Kernel>
  /*                                      */ Construct_Apollonius_vertex_2;
  typedef CGAL::Construct_Apollonius_site_2<Kernel>
  /*                                        */ Construct_Apollonius_site_2;


  // PREDICATES
  //-----------
  typedef CGAL::Ag2_compare_x_2<Kernel>                 Compare_x_2;
  typedef CGAL::Ag2_compare_y_2<Kernel>                 Compare_y_2;
  typedef CGAL::Ag2_compare_weight_2<Kernel>            Compare_weight_2;
  typedef CGAL::AG2_Orientation_test_2<Kernel>          Orientation_2;
  typedef CGAL::Is_hidden_2<Kernel,MTag>                Is_hidden_2;
  typedef CGAL::Oriented_side_of_bisector_2<Kernel,MTag> 
  /*                                          */ Oriented_side_of_bisector_2;
  typedef CGAL::Vertex_conflict_2<Kernel,MTag>          Vertex_conflict_2;
  typedef CGAL::Finite_edge_interior_conflict_2<Kernel,MTag>
  /*                                      */ Finite_edge_interior_conflict_2;
  typedef CGAL::Infinite_edge_interior_conflict_2<Kernel,MTag>
  /*                                    */ Infinite_edge_interior_conflict_2;
  typedef CGAL::Is_degenerate_edge_2<Kernel,MTag>       Is_degenerate_edge_2;


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

CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_TRAITS_2_H

// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Segment_Voronoi_diagram_traits_2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_2_H


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
//                           are same points
//-----------------------------------------------------------------------


template< class K >
class Are_same_points_2
{
public:
  typedef typename K::Point_2      Point_2;

  //  typedef typename K::Compare_x_2  compare_x_2;
  //  typedef typename K::Compare_y_2  compare_y_2;

  typedef bool                     result_type;

  struct Arity {};

public:

  bool operator()(const Point_2& p, const Point_2& q) const
  {
    return (p == q);
  }
};

//-----------------------------------------------------------------------
//                       oriented side of bisector
//-----------------------------------------------------------------------

template<class K, class Method_tag>
inline Comparison_result
svd_compare_distance_2(const typename K::Point_2& q,
		       const typename K::Segment_2& s,
		       const typename K::Point_2& p,
		       Cartesian_tag, Method_tag method_tag)
{
  return svd_compare_distance_ftC2(q.x(), q.y(),
				   s[0].x(), s[0].y(),
				   s[1].x(), s[1].y(),
				   p.x(), p.y(), method_tag);
}

template<class K, class Method_tag>
inline Comparison_result
svd_compare_distance_2(const typename K::Point_2& q,
		       const typename K::Segment_2& s,
		       const typename K::Point_2& p,
		       Homogeneous_tag, Method_tag method_tag)
{
  return svd_compare_distanceH2(q.hx(), q.hy(), q.hw(),
				s[0].hx(), s[0].hy(), s[0].hw(),
				s[1].hx(), s[1].hy(), s[1].hw(),
				p.hx(), p.hy(), p.hw(), method_tag);
}


template<class K, class Method_tag>
inline Comparison_result
svd_compare_distance_2(const typename K::Point_2& q,
		       const typename K::Segment_2& s1,
		       const typename K::Segment_2& s2,
		       Cartesian_tag, Method_tag method_tag)
{
  return svd_compare_distance_ftC2(q.x(), q.y(),
				   s1[0].x(), s1[0].y(),
				   s1[1].x(), s1[1].y(),
				   s2[0].x(), s2[0].y(),
				   s2[1].x(), s2[1].y(), method_tag);
}

template<class K, class Method_tag>
inline Comparison_result
svd_compare_distance_2(const typename K::Point_2& q,
		       const typename K::Segment_2& s1,
		       const typename K::Segment_2& s2,
		       Homogeneous_tag, Method_tag method_tag)
{
  return svd_compare_distanceH2(q.hx(), q.hy(), q.hw(),
				s1[0].hx(), s1[0].hy(), s1[0].hw(),
				s1[1].hx(), s1[1].hy(), s1[1].hw(),
				s2[0].hx(), s2[0].hy(), s2[0].hw(),
				s2[1].hx(), s2[1].hy(), s2[1].hw(),
				method_tag);
}


template<class K, class Method_tag>
class Svd_oriented_side_of_bisector_2
{
public:
  typedef typename K::Site_2     Site_2;
  typedef Oriented_side          result_type;

  typedef typename K::Point_2    Point_2;
  typedef typename K::Segment_2  Segment_2;
  typedef typename K::Rep_tag    Rep_tag;

  struct Arity {};

private:
  inline Comparison_result
  operator()(const Point_2& q,
	     const Point_2& p1, const Point_2& p2) const
  {
    CGAL_precondition( p1 != p2 );

    if ( q == p1 ) { return SMALLER; }
    if ( q == p2 ) { return LARGER; }
    
    return compare_distance_to_point(q, p1, p2);
  }

  inline Comparison_result
  operator()(const Point_2& q,
	     const Point_2& p, const Segment_2& s) const
  {
    CGAL_precondition( !s.is_degenerate() );

    return
      opposite( svd_compare_distance_2<K>(q, s, p, Rep_tag(), Method_tag()) );
  }

  inline Comparison_result
  operator()(const Point_2& q,
	     const Segment_2& s, const Point_2& p) const
  {
    if ( q == p ) { return LARGER; }
    if ( q == s.source() ) { return SMALLER; }
    if ( q == s.target() ) { return SMALLER; }

    return svd_compare_distance_2<K>(q, s, p, Rep_tag(), Method_tag());
  }

  inline Comparison_result
  operator()(const Point_2& q,
	     const Segment_2& s1, const Segment_2& s2) const
  {
    CGAL_precondition( !s1.is_degenerate() );
    CGAL_precondition( !s2.is_degenerate() );

    if (  ( q == s1.source() && q == s2.source() ) ||
	  ( q == s1.source() && q == s2.target() ) ||
	  ( q == s1.target() && q == s2.source() ) ||
	  ( q == s1.target() && q == s2.target() )  ) {
      return EQUAL;
    }

    if (  ( q == s1.source() || q == s1.target() ) &&
	  ( q != s2.source() && q != s2.target() )  ) {
      return SMALLER;
    }

    if (  ( q == s2.source() || q == s2.target() ) &&
	  ( q != s1.source() && q != s1.target() )  ) {
      return LARGER;
    }

    if ( (s1.source() == s2.source() && s1.target() == s2.target()) ||
	 (s1.source() == s2.target() && s1.target() == s2.source()) ) {
      return EQUAL;
    }

    return svd_compare_distance_2<K>(q, s1, s2, Rep_tag(), Method_tag());
  }

public:
  inline Oriented_side
  operator()(const Site_2& t1, const Site_2& t2,
	     const Point_2& q) const
  {
    Comparison_result r;

    if ( t1.is_point() && t2.is_point() ) {
      r = operator()(q, t1.point(), t2.point());
      if ( r == LARGER ) { return ON_NEGATIVE_SIDE; }
      if ( r == SMALLER ) { return ON_POSITIVE_SIDE; }
      return ON_ORIENTED_BOUNDARY;
    } else if ( t1.is_segment() && t2.is_point() ) {
      r = operator()(q, t1.segment(), t2.point());
      if ( r == LARGER ) { return ON_NEGATIVE_SIDE; }
      if ( r == SMALLER ) { return ON_POSITIVE_SIDE; }
      return ON_NEGATIVE_SIDE;
    } else if ( t1.is_point() && t2.is_segment() ) {
      r = operator()(q, t1.point(), t2.segment());
      if ( r == LARGER ) { return ON_NEGATIVE_SIDE; }
      if ( r == SMALLER ) { return ON_POSITIVE_SIDE; }
      return ON_POSITIVE_SIDE;
    } else {
      r = operator()(q, t1.segment(), t2.segment());
      if ( r == LARGER ) { return ON_NEGATIVE_SIDE; }
      if ( r == SMALLER ) { return ON_POSITIVE_SIDE; }
      return ON_ORIENTED_BOUNDARY;
    }
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
  Sign operator()(const Site_2& p, const Site_2& q,
		  const Site_2& t) const
  {
    return
      svd_vertex_conflict_ftC2<K,Method_tag>(q, p, t, Method_tag());
  }

  inline
  Sign operator()(const Site_2& p, const Site_2& q,
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
  
  inline
  bool operator()(const Site_2& p, const Site_2& q,
		  const Site_2& t, Sign sgn) const
  {
    return
      svd_finite_edge_conflict_ftC2<K,Method_tag>(p, q, t,
						  sgn, Method_tag());
  }


  inline
  bool operator()(const Site_2& p, const Site_2& q,
		  const Site_2& r, const Site_2& t, Sign sgn) const
  {
    return
      svd_finite_edge_conflict_ftC2<K,Method_tag>(p, q, r, t,
						  sgn, Method_tag());
  }

  inline
  bool operator()(const Site_2& p, const Site_2& q,
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

  struct Arity {};

public:
  
  inline
  bool operator()(const Site_2& q, const Site_2& r, const Site_2& s,
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

  struct Arity {};

public:
  
  inline
  bool operator()(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& s) const
  {
    return
      svd_is_degenerate_edge_ftC2<K,Method_tag>(p, q, r, s, Method_tag());
  }
};

//-----------------------------------------------------------------------
//                          do intersect
//-----------------------------------------------------------------------

template< class K, class Method_tag >
class Svd_do_intersect_2
{
public:
  typedef typename K::Segment_2  Segment_2;
  typedef std::pair<int,int>     result_type;

  struct Arity {};

public:
  
  result_type
  operator()(const Segment_2& s1, const Segment_2& s2) const
  {
    return
      svd_do_intersect_ftC2<K,Method_tag>(s1, s2, Method_tag());
  }
};




//-----------------------------------------------------------------------
// the Traits class
//-----------------------------------------------------------------------
template < class R, class MTag = Sqrt_field_tag >
class Segment_Voronoi_diagram_traits_2
{
public:
  //-----------------------------------------------------------------------
  //                  TYPE DEFINITIONS
  //-----------------------------------------------------------------------

  // BASIC TYPES
  //------------
  
  typedef Segment_Voronoi_diagram_kernel_wrapper_2<R>   Kernel;
  typedef Kernel                                        K;
  typedef R                                             Rep;
  typedef MTag                                          Method_tag;

  typedef typename Kernel::Point_2                Point_2;
  typedef typename Kernel::Line_2                 Line_2;
  typedef typename Kernel::Segment_2              Segment_2;
  typedef typename Kernel::Ray_2                  Ray_2;
  typedef typename Kernel::Circle_2               Circle_2;

  typedef Segment_Voronoi_diagram_site_2<R>       Site_2;

  typedef typename Kernel::FT                     FT;
  typedef typename Kernel::RT                     RT;

  typedef typename Kernel::Rep_tag                Rep_tag;

public:
  // CONSTRUCTIONS
  //--------------
  // vertex and Voronoi circle
  typedef CGAL::Construct_svd_vertex_2<K,MTag>  Construct_svd_vertex_2;

#if 0
  typedef CGAL::Construct_svd_circle_2<K,MTag>  Construct_svd_circle_2;
#endif

  // PREDICATES
  //-----------
  typedef typename Kernel::Compare_x_2                  Compare_x_2;
  typedef typename Kernel::Compare_y_2                  Compare_y_2;
  typedef typename Kernel::Orientation_2                Orientation_2;
  typedef CGAL::Are_same_points_2<K>                    Are_same_points_2;
  typedef CGAL::Svd_oriented_side_of_bisector_2<K,MTag> 
  /*                                          */ Oriented_side_of_bisector_2;
  typedef CGAL::Svd_vertex_conflict_2<K,MTag >             Vertex_conflict_2;
  typedef CGAL::Svd_finite_edge_interior_conflict_2<K,MTag >
  /*                                      */ Finite_edge_interior_conflict_2;
  typedef CGAL::Svd_infinite_edge_interior_conflict_2<K,MTag>
  /*                                    */ Infinite_edge_interior_conflict_2;
  typedef CGAL::Svd_is_degenerate_edge_2<K,MTag>        Is_degenerate_edge_2;
  typedef CGAL::Svd_do_intersect_2<K,MTag>              Do_intersect_2;


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

#if 0
  Construct_svd_circle_2
  construct_svd_circle_2_object() const {
    return Construct_svd_circle_2();
  }
#endif

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

  Do_intersect_2
  do_intersect_2_object() const {
    return Do_intersect_2();
  }

};


CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_2_H

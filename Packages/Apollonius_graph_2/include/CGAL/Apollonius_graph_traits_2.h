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
// file          : include/CGAL/Apollonius_graph_traits_2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_APOLLONIUS_GRAPH_TRAITS_2_H
#define CGAL_APOLLONIUS_GRAPH_TRAITS_2_H

#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Parabola_segment_2.h>
#include <CGAL/Hyperbola_2.h>
#include <CGAL/Hyperbola_segment_2.h>
#include <CGAL/Hyperbola_ray_2.h>

#include <CGAL/Filtered_kernel.h>
#include <CGAL/Filtered_predicate.h>

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/predicates/Apollonius_graph_ftC2.h>
#include <CGAL/Apollonius_graph_constructions_ftC2.h>
#endif

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/predicates/Apollonius_graph_rtH2.h>
#include <CGAL/Apollonius_graph_constructions_rtH2.h>
#endif


#include <CGAL/Kernel_traits.h>
#include <CGAL/Number_type_traits.h>

#include <CGAL/Apollonius_graph_kernel_wrapper_2.h>
#include <CGAL/Apollonius_graph_constructions_C2.h>


#define KEEP_MOST_TYPES_IN_TRAITS 0


CGAL_BEGIN_NAMESPACE


//***********************************************************************
//***********************************************************************
//                              PREDICATES
//***********************************************************************
//***********************************************************************

//-----------------------------------------------------------------------
//                        Compare weight
//-----------------------------------------------------------------------

template < class K >
class Compare_weight_2
{
public:
  typedef typename K::Site_2  Site_2;
  typedef Comparison_result             result_type;

  inline
  Comparison_result operator()(const Site_2& p,
			       const Site_2& q) const
  {
    return CGAL_NTS compare(p.weight(), q.weight());
  }
};

//-----------------------------------------------------------------------
//                        Is hidden
//-----------------------------------------------------------------------

template < class K >
inline
bool
ad_is_hidden_test_2(const typename K::Site_2& p,
		    const typename K::Site_2& q,
		    Cartesian_tag, Sqrt_field_tag )
{
  return ad_is_hidden_test_sqrtf_C2(p.x(), p.y(), p.weight(),
				    q.x(), q.y(), q.weight());
}


template < class K >
inline
bool
ad_is_hidden_test_2(const typename K::Site_2& p,
		    const typename K::Site_2& q,
		    Cartesian_tag, Ring_tag )
{
  return ad_is_hidden_test_ring_C2(p.x(), p.y(), p.weight(),
				   q.x(), q.y(), q.weight());
}



template < class K, class Method_tag >
inline
bool
ad_is_hidden_test_2(const typename K::Site_2& p,
		    const typename K::Site_2& q,
		    Homogeneous_tag)
{
  Sign s = sign_of_ad_distance2_testH2(p.hx(), p.hy(), p.hw(), 
				       p.weight(),
				       q.hx(), q.hy(), q.hw(),
				       q.weight());
  if ( s == POSITIVE ) { return false; }
  return (CGAL_NTS compare(p.weight(), q.weight()) != SMALLER);
}

template < class K, class Method_tag >
inline
bool
ad_is_hidden_test_2(const typename K::Site_2& p,
		    const typename K::Site_2& q,
		    Cartesian_tag tag)
{
  return ad_is_hidden_test_2< K >(p, q, tag, Method_tag());
}


template< class K, class Method_tag >
class Is_hidden_2
{
public:
  typedef typename K::Site_2   Site_2;
  typedef bool                           result_type;

  inline bool operator()(const Site_2 &p,
			 const Site_2 &q) const
  {
    typedef typename K::Rep_tag Tag;
    return ad_is_hidden_test_2<K,Method_tag>(p, q, Tag());
  }
};


//-----------------------------------------------------------------------
//                    Oriented side of bisector
//-----------------------------------------------------------------------


template < class K >
inline
Comparison_result
ad_distances_test_2(const typename K::Site_2& p1,
		    const typename K::Site_2& p2,
		    const typename K::Point_2& p,
		    Cartesian_tag, Sqrt_field_tag )
{
  return
    compare_ad_distances_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
				       p2.x(), p2.y(), p2.weight(),
				       p.x(),  p.y());
}


template < class K >
inline
Comparison_result
ad_distances_test_2(const typename K::Site_2& p1,
		    const typename K::Site_2& p2,
		    const typename K::Point_2& p,
		    Cartesian_tag, Ring_tag)
{
  return compare_ad_distances_test_ring_C2(p1.x(), p1.y(), p1.weight(),
					   p2.x(), p2.y(), p2.weight(),
					   p.x(),  p.y());
}



template < class K, class Method_tag >
inline
Comparison_result
ad_distances_test_2(const typename K::Site_2& p1,
		    const typename K::Site_2& p2,
		    const typename K::Point_2& p, Cartesian_tag tag)
{
  return ad_distances_test_2< K >(p1, p2, p, tag, Method_tag());
}



template < class K, class Method_tag >
inline
Comparison_result
ad_distances_test_2(const typename K::Site_2& p1,
		    const typename K::Site_2& p2,
		    const typename K::Point_2& p, Homogeneous_tag )
{
  return compare_ad_distances_testH2(p1.hx(), p1.hy(), p1.hw(),
				     p1.weight(),
				     p2.hx(), p2.hy(), p2.hw(),
				     p2.weight(),
				     p.hx(), p.hy(), p.hw());
}





template < class K, class Method_tag >
inline
Comparison_result
ad_distances_test_2(const typename K::Site_2& p1,
		    const typename K::Site_2& p2,
		    const typename K::Point_2& p)
{
  typedef typename K::Rep_tag Tag;
  return ad_distances_test_2<K,Method_tag>(p1, p2, p, Tag());
}



template< class K, class Method_tag >
class Oriented_side_of_bisector_2
{
public:
  typedef typename K::Point_2             Point_2;
  typedef typename K::Site_2    Site_2;
  typedef Oriented_side                   result_type;

  inline Oriented_side operator()(const Site_2& p1,
				  const Site_2& p2,
				  const Point_2 &p) const
  {
    Comparison_result r = ad_distances_test_2<K,Method_tag>(p1, p2, p);

    if ( r == EQUAL ) { return ON_ORIENTED_BOUNDARY; }
    return ( r == LARGER ) ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE;
  }
};



//-----------------------------------------------------------------------
//                        Vertex conflict
//-----------------------------------------------------------------------


template < class K >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2&  q,
		   Cartesian_tag, Sqrt_field_tag )
{
  return ad_incircle_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
				   p2.x(), p2.y(), p2.weight(),
				    q.x(),  q.y(),  q.weight());
}


template < class K >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2&  q,
		   Cartesian_tag, Ring_tag )
{
  return ad_incircle_test_ring_C2(p1.x(), p1.y(), p1.weight(),
				  p2.x(), p2.y(), p2.weight(),
				   q.x(),  q.y(),  q.weight());
}


template < class K, class Method_tag >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2&  q,
		   Cartesian_tag tag)
{
  return ad_incircle_test_2< K >(p1, p2, q, tag, Method_tag());
}


template < class K, class Method_tag >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2& q,
		   Homogeneous_tag )
{
  return 
    ad_incircle_testH2(p1.hx(), p1.hy(), p1.hw(), p1.weight(),
		       p2.hx(), p2.hy(), p2.hw(), p2.weight(),
		        q.hx(),  q.hy(),  q.hw(),  q.weight());
}


//-----------------------------------------------------------------------


template < class K >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2& p3,
		   const typename K::Site_2&  q,
		   Cartesian_tag, Sqrt_field_tag )
{
  return ad_incircle_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
				   p2.x(), p2.y(), p2.weight(),
				   p3.x(), p3.y(), p3.weight(),
				    q.x(),  q.y(),  q.weight());
}


template < class K >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2& p3,
		   const typename K::Site_2&  q,
		   Cartesian_tag, Ring_tag )
{
  return ad_incircle_test_ring_C2(p1.x(), p1.y(), p1.weight(),
				  p2.x(), p2.y(), p2.weight(),
				  p3.x(), p3.y(), p3.weight(),
				  q.x(),  q.y(),   q.weight());
}



template < class K, class Method_tag >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2& p3,
		   const typename K::Site_2&  q,
		   Cartesian_tag tag)
{
  return ad_incircle_test_2< K >(p1, p2, p3, q, tag, Method_tag());
}


template < class K, class Method_tag >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2& p3,
		   const typename K::Site_2& q,
		   Homogeneous_tag )
{
  return 
    ad_incircle_testH2(p1.hx(), p1.hy(), p1.hw(), p1.weight(),
		       p2.hx(), p2.hy(), p2.hw(), p2.weight(),
		       p3.hx(), p3.hy(), p3.hw(), p3.weight(),
		        q.hx(),  q.hy(),  q.hw(),  q.weight());
}


template < class K, class Method_tag >
class Vertex_conflict_2
{
public:
  typedef typename K::Site_2      Site_2;
  typedef Sign                               result_type;

  inline
  Sign operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& q) const
  {
    typedef typename K::Rep_tag Tag;
    return ad_incircle_test_2<K,Method_tag>(p1, p2, p3, q, Tag());
  }


  inline
  Sign operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& q) const
  {
    typedef typename K::Rep_tag Tag;
    return ad_incircle_test_2<K,Method_tag>(p1, p2, q, Tag());
  }
 


};

//-----------------------------------------------------------------------
//                    Finite edge interior conflict
//-----------------------------------------------------------------------

template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Sqrt_field_tag)
{
  return
    ad_finite_edge_test_degenerated_sqrtf_C2(p1.x(), p1.y(),
					     p1.weight(),
					     p2.x(), p2.y(),
					     p2.weight(),
					     q.x(),  q.y(),
					     q.weight(), b);
}

template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Ring_tag)
{
  return
    ad_finite_edge_test_degenerated_ring_C2(p1.x(), p1.y(),
					    p1.weight(),
					    p2.x(), p2.y(),
					    p2.weight(),
					    q.x(),  q.y(),
					    q.weight(), b);
}


template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag tag)
{
  return
    ad_finite_edge_test_2< K >(p1, p2, q, b, tag, Method_tag());
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& q,
		      bool b, Homogeneous_tag)
{
  return
    ad_finite_edge_test_degeneratedH2(p1.hx(), p1.hy(),
				      p1.hw(),
				      p1.weight(),
				      p2.hx(), p2.hy(),
				      p2.hw(),
				      p2.weight(),
				      q.hx(),  q.hy(),
				      q.hw(),
				      q.weight(), b);
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& q, bool b)
{
  typedef typename K::Rep_tag Tag;
  return ad_finite_edge_test_2<K,Method_tag>(p1, p2, q, b, Tag());
}

//-----------------------------------------------------------------------

template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Sqrt_field_tag)
{
  return ad_finite_edge_test_degenerated_sqrtf_C2(p1.x(), p1.y(),
						  p1.weight(),
						  p2.x(), p2.y(),
						  p2.weight(),
						  p3.x(), p3.y(),
						  p3.weight(),
						  q.x(),  q.y(),
						  q.weight(), b);
}

template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Ring_tag)
{
  return ad_finite_edge_test_degenerated_ring_C2(p1.x(), p1.y(),
						 p1.weight(),
						 p2.x(), p2.y(),
						 p2.weight(),
						 p3.x(), p3.y(),
						 p3.weight(),
						 q.x(),  q.y(),
						 q.weight(), b);
}


template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag tag)
{
  return
    ad_finite_edge_test_2< K >(p1, p2, p3, q, b, tag, Method_tag());
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& q,
		      bool b, Homogeneous_tag)
{
  return
    ad_finite_edge_test_degeneratedH2(p1.hx(), p1.hy(),
				      p1.hw(),
				      p1.weight(),
				      p2.hx(), p2.hy(),
				      p2.hw(),
				      p2.weight(),
				      p3.hx(), p3.hy(),
				      p3.hw(),
				      p3.weight(),
				      q.hx(),  q.hy(),
				      q.hw(),
				      q.weight(), b);
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& q, bool b)
{
  typedef typename K::Rep_tag Tag;
  return ad_finite_edge_test_2<K,Method_tag>(p1, p2, p3, q, b, Tag());
}

//-----------------------------------------------------------------------


template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& p4,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Sqrt_field_tag)
{
  return
    ad_finite_edge_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
				 p2.x(), p2.y(), p2.weight(),
				 p3.x(), p3.y(), p3.weight(),
				 p4.x(), p4.y(), p4.weight(),
				 q.x(),  q.y(), q.weight(), b);
}

template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& p4,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Ring_tag)
{
  return
    ad_finite_edge_test_ring_C2(p1.x(), p1.y(), p1.weight(),
				p2.x(), p2.y(), p2.weight(),
				p3.x(), p3.y(), p3.weight(),
				p4.x(), p4.y(), p4.weight(),
				 q.x(),  q.y(),  q.weight(), b);
}


template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& p4,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag tag)
{
  return ad_finite_edge_test_2< K >(p1, p2, p3, p4, q, b,
				    tag, Method_tag());
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& p4,
		      const typename K::Site_2& q,
		      bool b, Homogeneous_tag)
{
  return
    ad_Voronoi_diagram_finite_edge_testH2(p1.hx(), p1.hy(), p1.hw(),
					  p1.weight(),
					  p2.hx(), p2.hy(), p2.hw(),
					  p2.weight(),
					  p3.hx(), p3.hy(), p3.hw(),
					  p3.weight(),
					  p4.hx(), p4.hy(), p4.hw(),
					  p4.weight(),
					  q.hx(),  q.hy(), q.hw(),
					  q.weight(), b);
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& p4,
		      const typename K::Site_2& q,
		      bool b)
{
  typedef typename K::Rep_tag Tag;
  return ad_finite_edge_test_2<K,Method_tag>
    (p1, p2, p3, p4, q, b, Tag());
}




template < class K, class Method_tag >
class Finite_edge_interior_conflict_2
{
public:
  typedef typename K::Site_2  Site_2;
  typedef bool                           result_type;

  inline
  bool operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& q, bool b) const
  {
    return ad_finite_edge_test_2<K,Method_tag>(p1, p2, p3, q, b);
  }

  bool operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& q, bool b) const
  {
    return ad_finite_edge_test_2<K,Method_tag>(p1, p2, q, b);
  }



  inline
  bool operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& p4,
		  const Site_2& q,
		  bool b) const
  {
    return ad_finite_edge_test_2<K,Method_tag>(p1, p2, p3, p4, q, b);
  }
};


//-----------------------------------------------------------------------
//                   Infinite edge interior conflict
//-----------------------------------------------------------------------

template < class K >
inline
bool
ad_infinite_edge_test_2(const typename K::Site_2& p2,
			const typename K::Site_2& p3,
			const typename K::Site_2& p4,
			const typename K::Site_2& q,
			bool b, Cartesian_tag, Sqrt_field_tag)
{
  return
    ad_infinite_edge_test_sqrtf_C2(p2.x(), p2.y(), p2.weight(),
				   p3.x(), p3.y(), p3.weight(),
				   p4.x(), p4.y(), p4.weight(),
				    q.x(),  q.y(),  q.weight(), b);
}


template < class K >
inline
bool
ad_infinite_edge_test_2(const typename K::Site_2& p2,
			const typename K::Site_2& p3,
			const typename K::Site_2& p4,
			const typename K::Site_2& q,
			bool b, Cartesian_tag, Ring_tag)
{
  return
    ad_infinite_edge_test_ring_C2(p2.x(), p2.y(), p2.weight(),
				  p3.x(), p3.y(), p3.weight(),
				  p4.x(), p4.y(), p4.weight(),
				   q.x(),  q.y(),  q.weight(), b);
}


template < class K, class Method_tag >
inline
bool
ad_infinite_edge_test_2(const typename K::Site_2& p2,
			const typename K::Site_2& p3,
			const typename K::Site_2& p4,
			const typename K::Site_2& q,
			bool b, Cartesian_tag tag)
{
  return
    ad_infinite_edge_test_2<K>(p2, p3, p4, q, b, tag, Method_tag());
}



template < class K, class Method_tag >
inline
bool
ad_infinite_edge_test_2(const typename K::Site_2& p2,
			const typename K::Site_2& p3,
			const typename K::Site_2& p4,
			const typename K::Site_2& q,
			bool b, Homogeneous_tag)
{
  return
    ad_infinite_edge_testH2(p2.hx(), p2.hy(), p2.hw(),
			    p2.weight(),
			    p3.hx(), p3.hy(), p3.hw(),
			    p3.weight(),
			    p4.hx(), p4.hy(), p4.hw(),
			    p4.weight(),
			    q.hx(),  q.hy(), q.hw(),
			    q.weight(), b);
}

template < class K, class Method_tag >
inline
bool
ad_infinite_edge_test_2(const typename K::Site_2& p2,
			const typename K::Site_2& p3,
			const typename K::Site_2& p4,
			const typename K::Site_2& q, bool b)
{
  typedef typename K::Rep_tag Tag;
  return ad_infinite_edge_test_2<K,Method_tag>(p2, p3, p4, q, b, Tag());
}


template < class K, class Method_tag >
class Infinite_edge_interior_conflict_2
{
public:
  typedef typename K::Site_2 Site_2;
  typedef bool                          result_type;

  inline
  bool operator()(const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& p4,
		  const Site_2& q, bool b) const
  {
    return ad_infinite_edge_test_2<K,Method_tag>(p2, p3, p4, q, b);
  }
};


//-----------------------------------------------------------------------
//                          Is degenerate
//-----------------------------------------------------------------------



template < class K >
inline
bool
ad_is_degenerate_edge_test_2(const typename K::Site_2& p1,
			     const typename K::Site_2& p2,
			     const typename K::Site_2& p3,
			     const typename K::Site_2& p4,
			     Cartesian_tag, Sqrt_field_tag)
{
  return
    ad_is_degenerate_edge_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
					p2.x(), p2.y(), p2.weight(),
					p3.x(), p3.y(), p3.weight(),
					p4.x(), p4.y(), p4.weight());
}

template < class K >
inline
bool
ad_is_degenerate_edge_test_2(const typename K::Site_2& p1,
			     const typename K::Site_2& p2,
			     const typename K::Site_2& p3,
			     const typename K::Site_2& p4,
			     Cartesian_tag, Ring_tag)
{
  return
    ad_is_degenerate_edge_test_ring_C2(p1.x(), p1.y(), p1.weight(),
				       p2.x(), p2.y(), p2.weight(),
				       p3.x(), p3.y(), p3.weight(),
				       p4.x(), p4.y(), p4.weight());
}


template < class K, class Method_tag >
inline
bool
ad_is_degenerate_edge_test_2(const typename K::Site_2& p1,
			     const typename K::Site_2& p2,
			     const typename K::Site_2& p3,
			     const typename K::Site_2& p4,
			     Cartesian_tag tag)
{
  return
    ad_is_degenerate_edge_test_2< K >(p1, p2, p3, p4, tag, Method_tag());
}


template < class K, class Method_tag >
inline
bool
ad_is_degenerate_edge_test_2(const typename K::Site_2& p1,
			     const typename K::Site_2& p2,
			     const typename K::Site_2& p3,
			     const typename K::Site_2& p4,
			     Homogeneous_tag)
{
  return
    ad_is_degenerate_edge_testH2(p1.hx(), p1.hy(), p1.hw(),
				 p1.weight(),
				 p2.hx(), p2.hy(), p2.hw(),
				 p2.weight(),
				 p3.hx(), p3.hy(), p3.hw(),
				 p3.weight(),
				 p4.hx(), p4.hy(), p4.hw(),
				 p4.weight());
}

template < class K, class Method_tag >
inline
bool
ad_is_degenerate_edge_test_2(const typename K::Site_2& p1,
			     const typename K::Site_2& p2,
			     const typename K::Site_2& p3,
			     const typename K::Site_2& p4)
{
  typedef typename K::Rep_tag Tag;
  return
    ad_is_degenerate_edge_test_2<K,Method_tag>(p1, p2, p3, p4, Tag());
}



template < class K, class Method_tag >
class Is_degenerate_edge_2
{
public:
  typedef typename K::Site_2    Site_2;
  typedef bool                             result_type;

  inline
  bool operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& p4) const
  {
    return
      ad_is_degenerate_edge_test_2<K,Method_tag>(p1, p2, p3, p4);
  }
};




//-----------------------------------------------------------------------
// the Traits class
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
  typedef Apollonius_graph_kernel_wrapper_2<Rep>        Kernel;
  typedef Apollonius_graph_traits_2<Rep,MTag>           Self;

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

private:
  typedef typename Site_2::FT                           Weight;

public:
#if KEEP_MOST_TYPES_IN_TRAITS
  typedef CGAL::Parabola_segment_2<Point_2,Weight,Line_2>
  /*                                                 */ Parabola_segment_2;
  typedef CGAL::Hyperbola_2<Self>                       Hyperbola_2;
  typedef CGAL::Hyperbola_ray_2<Self>                   Hyperbola_ray_2;
  typedef CGAL::Hyperbola_segment_2<Self>
  /*                                                 */ Hyperbola_segment_2;
#endif

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


#if KEEP_MOST_TYPES_IN_TRAITS
  // bisectors and subsets
  typedef CGAL::Construct_Apollonius_bisector_2<Self>
  /*                                    */ Construct_Apollonius_bisector_2;
  typedef CGAL::Construct_Apollonius_bisector_ray_2<Self>
  /*                                */ Construct_Apollonius_bisector_ray_2;
  typedef CGAL::Construct_Apollonius_bisector_segment_2<Self>
  /*                            */ Construct_Apollonius_bisector_segment_2;

  // primal edges
  typedef CGAL::Construct_Apollonius_primal_ray_2<Self> 
  /*                                   */ Construct_Apollonius_primal_ray_2;
  typedef CGAL::Construct_Apollonius_primal_segment_2<Self>
  /*                               */ Construct_Apollonius_primal_segment_2;
#endif

  // PREDICATES
  //-----------
  typedef typename Kernel::Compare_x_2                  Compare_x_2;
  typedef typename Kernel::Compare_y_2                  Compare_y_2;
  typedef CGAL::Compare_weight_2<Kernel>                Compare_weight_2;
  typedef typename Kernel::Orientation_2                Orientation_2;
  typedef CGAL::Is_hidden_2<Kernel,MTag>                Is_hidden_2;
  typedef CGAL::Oriented_side_of_bisector_2<Kernel,MTag> 
  /*                                          */ Oriented_side_of_bisector_2;
  typedef CGAL::Vertex_conflict_2<Kernel,MTag >            Vertex_conflict_2;
  typedef CGAL::Finite_edge_interior_conflict_2<Kernel,MTag >
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

#if KEEP_MOST_TYPES_IN_TRAITS
  Construct_Apollonius_bisector_2
  construct_Apollonius_bisector_2_object() const {
    return Construct_Apollonius_bisector_2();
  }

  Construct_Apollonius_bisector_ray_2
  construct_Apollonius_bisector_ray_2_object() const {
    return Construct_Apollonius_bisector_ray_2();
  }

  Construct_Apollonius_bisector_segment_2
  construct_Apollonius_bisector_segment_2_object() const { 
    return Construct_Apollonius_bisector_segment_2(); 
  }

  Construct_Apollonius_primal_ray_2
  construct_Apollonius_primal_ray_2_object() const {
    return Construct_Apollonius_primal_ray_2(); 
  }

  Construct_Apollonius_primal_segment_2
  construct_Apollonius_primal_segment_2_object() const { 
    return Construct_Apollonius_primal_segment_2();
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


#if 0
// I AM REMOVING THE FILTERED KERNEL STUFF SO THAT I DON'T GET ANY
// PROBLEMS WITH PARTIAL SPECIALIZATION AND ALSO BECAUSE THE CURRENT
// DEFINITION OF FILTERED KERNEL DOES NOT SUPPORT MORE THAN ONE
// TEMPLATE ARGUMENT

// CHANGE THE TEMPLATE NAMES SO THAT THERE IS NO UNDERSCROEE IN THE
// BEGINNING

//-----------------------------------------------------------------------
// the Traits class for a filtered kernel
//-----------------------------------------------------------------------
template<class C_Traits,
	 class E_Traits =
	 Apollonius_graph_traits_2<Cartesian<MP_Float> > >
	 
class Apollonius_graph_filtered_traits_2
{
private:
  typedef typename C_Traits::R                             CK;
  typedef typename E_Traits::R                             EK;
  typedef typename Simple_cartesian<Interval_nt_advanced>  FK;

  typedef Apollonius_graph_traits_2<

  typedef Apollonius_graph_kernel_wrapper_2<_CK>  CK;
  typedef Apollonius_graph_kernel_wrapper_2<_EK>  EK;
  typedef Apollonius_graph_kernel_wrapper_2<_FK>  FK;
  typedef Extended_cartesian_converter<_C2E>  C2E;
  typedef Extended_cartesian_converter<_C2F>  C2F;

public:
  //-----------------------------------------------------------------------
  //                  TYPE DEFINITIONS
  //-----------------------------------------------------------------------

  // BASIC TYPES
  //------------
  typedef Kernel_wrapper_2< Filtered_kernel<CK,EK,FK,C2E,C2F> >  Kernel;
  typedef Kernel                                         Rep;
  typedef MTag                                           Method_tag;
  typedef typename Kernel::RT                            Weight;
  typedef typename Kernel::Point_2                       Point_2;
  typedef typename Kernel::Site_2                        Site_2;
  //
  typedef typename Kernel::Object_2                      Object_2;
  typedef typename Kernel::Line_2                        Line_2;
  typedef typename Kernel::Ray_2                         Ray_2;
  typedef typename Kernel::Segment_2                     Segment_2;
  typedef CGAL::Parabola_segment_2<Point_2,Weight,Line_2>
  /*                                                  */ Parabola_segment_2;
  typedef CGAL::Hyperbola_2<Point_2,Weight>         Hyperbola_2;
  typedef CGAL::Hyperbola_ray_2<Point_2,Weight>     Hyperbola_ray_2;
  typedef CGAL::Hyperbola_segment_2<Point_2,Weight> Hyperbola_segment_2;

public:
  // CONSTRUCTIONS
  //--------------
  // vertex and dual site
  typedef CGAL::Construct_Apollonius_vertex_2<Kernel>
  /*                                      */ Construct_Apollonius_vertex_2;
  typedef CGAL::Construct_Apollonius_site_2<Kernel>
  /*                                        */ Construct_Apollonius_site_2;

  // bisectors and subsets
  typedef CGAL::Construct_Apollonius_bisector_2<Kernel>
  /*                                    */ Construct_Apollonius_bisector_2;
  typedef CGAL::Construct_Apollonius_bisector_ray_2<Kernel>
  /*                                */ Construct_Apollonius_bisector_ray_2;
  typedef CGAL::Construct_Apollonius_bisector_segment_2<Kernel>
  /*                            */ Construct_Apollonius_bisector_segment_2;

  // primal edges
  typedef CGAL::Construct_Apollonius_primal_ray_2<Kernel>
  /*                                   */ Construct_Apollonius_primal_ray_2;
  typedef CGAL::Construct_Apollonius_primal_segment_2<Kernel>
  /*                               */ Construct_Apollonius_primal_segment_2;


private:
  // Predicates for the construction kernel
  typedef CGAL::Compare_weight_2<CK>                   CK_compare_weight_2;
  typedef CGAL::Is_hidden_2<CK,MTag>                   CK_is_hidden_2;
  typedef CGAL::Oriented_side_of_bisector_2<CK,MTag> 
  /*                                      */ CK_oriented_side_of_bisector_2;
  typedef CGAL::Vertex_conflict_2<CK,MTag >            CK_vertex_conflict_2;
  typedef CGAL::Finite_edge_interior_conflict_2<CK,MTag >
  /*                                  */ CK_finite_edge_interior_conflict_2;
  typedef CGAL::Infinite_edge_interior_conflict_2<CK,MTag>
  /*                                */ CK_infinite_edge_interior_conflict_2;
  typedef CGAL::Is_degenerate_edge_2<CK,MTag>       CK_is_degenerate_edge_2;

  // Predicates for the exact kernel
  typedef CGAL::Compare_weight_2<EK>                   EK_compare_weight_2;
  typedef CGAL::Is_hidden_2<EK,MTag>                   EK_is_hidden_2;
  typedef CGAL::Oriented_side_of_bisector_2<EK,MTag> 
  /*                                      */ EK_oriented_side_of_bisector_2;
  typedef CGAL::Vertex_conflict_2<EK,MTag >            EK_vertex_conflict_2;
  typedef CGAL::Finite_edge_interior_conflict_2<EK,MTag >
  /*                                  */ EK_finite_edge_interior_conflict_2;
  typedef CGAL::Infinite_edge_interior_conflict_2<EK,MTag>
  /*                                */ EK_infinite_edge_interior_conflict_2;
  typedef CGAL::Is_degenerate_edge_2<EK,MTag>       EK_is_degenerate_edge_2;

  // Predicates for the filtering kernel
  typedef CGAL::Compare_weight_2<FK>                   FK_compare_weight_2;
  typedef CGAL::Is_hidden_2<FK,MTag>                   FK_is_hidden_2;
  typedef CGAL::Oriented_side_of_bisector_2<FK,MTag> 
  /*                                      */ FK_oriented_side_of_bisector_2;
  typedef CGAL::Vertex_conflict_2<FK,MTag >            FK_vertex_conflict_2;
  typedef CGAL::Finite_edge_interior_conflict_2<FK,MTag >
  /*                                  */ FK_finite_edge_interior_conflict_2;
  typedef CGAL::Infinite_edge_interior_conflict_2<FK,MTag>
  /*                                */ FK_infinite_edge_interior_conflict_2;
  typedef CGAL::Is_degenerate_edge_2<FK,MTag>       FK_is_degenerate_edge_2;

public:
  // PREDICATES
  //-----------
  typedef typename Kernel::Compare_x_2                   Compare_x_2;
  typedef typename Kernel::Compare_y_2                   Compare_y_2;
  typedef Filtered_predicate<EK_compare_weight_2,
    FK_compare_weight_2,C2E,C2F> Compare_weight_2;
  typedef typename Kernel::Orientation_2                 Orientation_2;
  typedef Filtered_predicate<EK_is_hidden_2,
    FK_is_hidden_2,C2E,C2F> Is_hidden_2;
  typedef Filtered_predicate<EK_oriented_side_of_bisector_2,
    FK_oriented_side_of_bisector_2,C2E,C2F> Oriented_side_of_bisector_2;
  typedef Filtered_predicate<EK_vertex_conflict_2,
    FK_vertex_conflict_2,C2E,C2F> Vertex_conflict_2;
  typedef Filtered_predicate<EK_finite_edge_interior_conflict_2,
    FK_finite_edge_interior_conflict_2,C2E,C2F>
  Finite_edge_interior_conflict_2;
  typedef Filtered_predicate<EK_infinite_edge_interior_conflict_2,
    FK_infinite_edge_interior_conflict_2,C2E,C2F>
  Infinite_edge_interior_conflict_2;
  typedef Filtered_predicate<EK_is_degenerate_edge_2,
    FK_is_degenerate_edge_2,C2E,C2F> Is_degenerate_edge_2;

public:
  //-----------------------------------------------------------------------
  //                  ACCESS TO OBJECTS
  //-----------------------------------------------------------------------

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

  Construct_Apollonius_bisector_2
  construct_Apollonius_bisector_2_object() const {
    return Construct_Apollonius_bisector_2();
  }

  Construct_Apollonius_bisector_ray_2
  construct_Apollonius_bisector_ray_2_object() const {
    return Construct_Apollonius_bisector_ray_2();
  }

  Construct_Apollonius_bisector_segment_2
  construct_Apollonius_bisector_segment_2_object() const { 
    return Construct_Apollonius_bisector_segment_2(); 
  }

  Construct_Apollonius_primal_ray_2
  construct_Apollonius_primal_ray_2_object() const {
    return Construct_Apollonius_primal_ray_2(); 
  }

  Construct_Apollonius_primal_segment_2
  construct_Apollonius_primal_segment_2_object() const { 
    return Construct_Apollonius_primal_segment_2();
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
#endif


CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_TRAITS_2_H

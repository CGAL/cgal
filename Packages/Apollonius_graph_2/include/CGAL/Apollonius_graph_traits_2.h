// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/predicates/Apollonius_graph_rtH2.h>
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

template<class K>
inline
Orientation
ag2_orientation_test_2(const typename K::Site_2& p,
		       const typename K::Site_2& q,
		       const typename K::Site_2& r,
		       Cartesian_tag)
{
  return ag2_orientation_test_C2(p.x(), p.y(), p.weight(),
				 q.x(), q.y(), q.weight(),
				 r.x(), r.y(), r.weight());
}

template<class K>
inline
Orientation
ag2_orientation_test_2(const typename K::Site_2& p,
		       const typename K::Site_2& q,
		       const typename K::Site_2& r,
		       Homogeneous_tag)
{
  CGAL_assertion( false );
  return COLLINEAR;
}

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
    typedef typename K::Rep_tag Tag;
    return ag2_orientation_test_2<K>(p, q, r, Tag());
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
  return (CGAL::compare(p.weight(), q.weight()) != SMALLER);
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
  typedef bool                 result_type;
  typedef Arity_tag<2>         Arity;

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
  typedef typename K::Site_2              Site_2;
  typedef Oriented_side                   result_type;
  typedef Arity_tag<3>                    Arity;

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
  typedef Sign                    result_type;
  struct Arity {};

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
  typedef bool                result_type;
  struct Arity {};

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
  typedef typename K::Site_2   Site_2;
  typedef bool                 result_type;
  typedef Arity_tag<5>         Arity;

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
  typedef bool                  result_type;
  typedef Arity_tag<4>          Arity;

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

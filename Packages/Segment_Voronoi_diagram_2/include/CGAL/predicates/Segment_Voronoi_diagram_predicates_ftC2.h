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




#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_PREDICATES_FTC2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_PREDICATES_FTC2_H

#include <CGAL/determinant.h>
#include <CGAL/predicates/Segment_Voronoi_diagram_predicates_C2.h>
#include <CGAL/predicates/check_filter.h>
#include <CGAL/Segment_Voronoi_diagram_kernel_wrapper_2.h>


CGAL_BEGIN_NAMESPACE

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

template<class K>
void svd_predicate_push_back_C2(const typename K::Site_2& t,
				typename K::FT v[], unsigned int& k,
				char site_types[], unsigned int& j)
{
  unsigned int step(0);

  if ( t.is_point() ) {
    site_types[j] = 'p';
    if ( t.is_exact() ) {
      site_types[j+1] = 'e';
      v[k] = t.point().x();
      v[k+1] = t.point().y();
      step = 2;
    } else {
      site_types[j+1] = 'i';
      typename K::Segment_2 s1 = t.supporting_segment(0);
      typename K::Segment_2 s2 = t.supporting_segment(1);
      v[k] = s1.source().x();
      v[k+1] = s1.source().y();
      v[k+2] = s1.target().x();
      v[k+3] = s1.target().y();
      v[k+4] = s2.source().x();
      v[k+5] = s2.source().y();
      v[k+6] = s2.target().x();
      v[k+7] = s2.target().y();
      step = 8;
    }
  } else {
    site_types[j] = 's';
    if ( t.is_exact() ) {
      site_types[j+1] = 'e';
      v[k] = t.source().x();
      v[k+1] = t.source().y();
      v[k+2] = t.target().x();
      v[k+3] = t.target().y();
      step = 4;
    } else {
      typename K::Segment_2 supp = t.supporting_segment();
      v[k] = supp.source().x();
      v[k+1] = supp.source().y();
      v[k+2] = supp.target().x();
      v[k+3] = supp.target().y();

      typename K::Segment_2 cs, cs2;
      char stype;
      if ( t.is_exact(0) ) {
	stype = '0';
	cs = t.crossing_segment(1);
      } else if ( t.is_exact(1) ) {
	stype = '1';
	cs = t.crossing_segment(0);
      } else {
	stype = 'i';
	cs = t.crossing_segment(0);
	cs2 = t.crossing_segment(1);
      }

      site_types[j+1] = stype;
      v[k+4] = cs.source().x();
      v[k+5] = cs.source().y();
      v[k+6] = cs.target().x();
      v[k+7] = cs.target().y();

      step = 8;

      if ( stype == 'i' ) {
	v[k+8] = cs2.source().x();
	v[k+9] = cs2.source().y();
	v[k+10] = cs2.target().x();
	v[k+11] = cs2.target().y();
	step = 12;
      }
    }
  }

  j += 2;
  k += step;
}


template<class K>
typename K::Site_2
get_site(const typename K::FT v[], unsigned int& k,
	 const char site_types[], unsigned int& j)
{
  typedef typename K::Point_2             Point_2;
  typedef typename K::Segment_2           Segment_2;
  typedef typename K::Site_2              Site_2;

  Site_2 t;

  unsigned int step(0);

  if ( site_types[j] == 'p' ) {
    if ( site_types[j+1] == 'e' ) {
      Point_2 p(v[k], v[k+1]);
      t = Site_2::construct_site_2(p);
      step = 2;
    } else {
      Point_2 p1(v[k], v[k+1]), p2(v[k+2], v[k+3]);
      Point_2 p3(v[k+4], v[k+5]), p4(v[k+6], v[k+7]);
      t = Site_2::construct_site_2(p1, p2, p3, p4);
      step = 8;
    }
  } else {
    if ( site_types[j+1] == 'e' ) {
      Point_2 p1(v[k], v[k+1]), p2(v[k+2], v[k+3]);
      t = Site_2::construct_site_2(p1, p2);
      step = 4;
    } else {
      if ( site_types[j+1] != 'i' ) {
	Point_2 p1(v[k], v[k+1]), p2(v[k+2], v[k+3]);
	Point_2 p3(v[k+4], v[k+5]), p4(v[k+6], v[k+7]);
	t = Site_2::construct_site_2(p1, p2, p3, p4, (site_types[j+1] == '0'));
	step = 8;
      } else {
	Point_2 p1(v[k], v[k+1]), p2(v[k+2], v[k+3]);
	Point_2 p3(v[k+4], v[k+5]), p4(v[k+6], v[k+7]);
	Point_2 p5(v[k+8], v[k+9]), p6(v[k+10], v[k+11]);
	t = Site_2::construct_site_2(p1, p2, p3, p4, p5, p6);
	step = 12;
      }
    }
  }

  j += 2;
  k += step;

  return t;
}

//--------------------------------------------------------------------------

template<typename Result_t, class Predicate, unsigned int Arity>
struct Svd_predicate_caller;

template<typename Result_t, class Predicate>
struct Svd_predicate_caller<Result_t, Predicate, 2>
{
  template<class S>
  Result_t operator()(const S t[]) const
    {
      return Predicate()(t[0], t[1]);
    }
};

template<typename Result_t, class Predicate>
struct Svd_predicate_caller<Result_t, Predicate, 3>
{
  template<class S>
  Result_t operator()(const S t[]) const
    {
      return Predicate()(t[0], t[1], t[2]);
    }

  template<class S, typename Data>
  Result_t operator()(const S t[], Data data) const
    {
      return Predicate()(t[0], t[1], t[2], data);
    }
};

template<typename Result_t, class Predicate>
struct Svd_predicate_caller<Result_t, Predicate, 4>
{
  template<class S>
  Result_t operator()(const S t[]) const
    {
      return Predicate()(t[0], t[1], t[2], t[3]);
    }

  template<class S, typename Data>
  Result_t operator()(const S t[], Data data) const
    {
      return Predicate()(t[0], t[1], t[2], t[3], data);
    }
};

template<typename Result_t, class Predicate>
struct Svd_predicate_caller<Result_t, Predicate, 5>
{
  template<class S>
  Result_t operator()(const S t[]) const
    {
      return Predicate()(t[0], t[1], t[2], t[3], t[4]);
    }

  template<class S, typename Data>
  Result_t operator()(const S t[], Data data) const
    {
      return Predicate()(t[0], t[1], t[2], t[3], t[4], data);
    }
};


//--------------------------------------------------------------------------

template<template<class Kernel> class Predicate_t,
	 typename Return_t, class FT, class ITag,
	 unsigned int Num_sites>
Return_t
svd_predicate_ftC2(const FT v[], const char site_types[])
{
  typedef Simple_cartesian<FT>                 Rep;
  typedef
  CGAL::Segment_Voronoi_diagram_kernel_wrapper_2<Rep,ITag>  Kernel;

  typedef typename Kernel::Site_2                   Site_2;
  typedef Predicate_t<Kernel>                       Predicate;

  typedef Svd_predicate_caller<Return_t, Predicate, Num_sites> Caller;

  must_be_filtered(FT());


  Site_2 t[Num_sites];

  for (unsigned int i = 0, k = 0, j = 0; i < Num_sites; i++) {
    t[i] = get_site<Kernel>(v, k, site_types, j);
  }

  Return_t result = Caller()(t);

  return result;
}

template<template<class Kernel, class MTag> class Predicate_t,
	 typename Return_t, class FT, class Method_tag,
	 class ITag, unsigned int Num_sites>
Return_t
svd_predicate_ftC2(const FT v[], const char site_types[])
{
  typedef Simple_cartesian<FT>                 Rep;
  typedef
  CGAL::Segment_Voronoi_diagram_kernel_wrapper_2<Rep,ITag>  Kernel;

  typedef typename Kernel::Site_2                   Site_2;
  typedef Predicate_t<Kernel,Method_tag>            Predicate;

  typedef Svd_predicate_caller<Return_t, Predicate, Num_sites> Caller;

  must_be_filtered(FT());


  Site_2 t[Num_sites];

  for (unsigned int i = 0, k = 0, j = 0; i < Num_sites; i++) {
    t[i] = get_site<Kernel>(v, k, site_types, j);
  }

  Return_t result = Caller()(t);

  return result;
}

template<template<class Kernel, class MTag> class Predicate_t,
	 typename Return_t, class FT, class Method_tag, 
	 class ITag, typename Data, unsigned int Num_sites>
Return_t
svd_predicate_ftC2(const FT v[], const char site_types[], Data data)
{
  typedef Simple_cartesian<FT>                 Rep;
  typedef
  CGAL::Segment_Voronoi_diagram_kernel_wrapper_2<Rep,ITag>  Kernel;

  typedef typename Kernel::Site_2                   Site_2;
  typedef Predicate_t<Kernel,Method_tag>            Predicate;

  typedef Svd_predicate_caller<Return_t, Predicate, Num_sites> Caller;

  must_be_filtered(FT());


  Site_2 t[Num_sites];

  for (unsigned int i = 0, k = 0, j = 0; i < Num_sites; i++) {
    t[i] = get_site<Kernel>(v, k, site_types, j);
  }

  Return_t result = Caller()(t, data);

  return result;
}

//--------------------------------------------------------------------------

template<template<class Kernel> class Predicate,
	 typename Return_t, class K,
	 unsigned int Num_sites>
Return_t
svd_predicate_C2(const typename K::Site_2 t[])
{
  typedef typename K::FT                 FT;
  typedef typename K::Intersections_tag  ITag;

  FT v[Num_sites * 12];
  char site_types[Num_sites * 2];

  for (unsigned int i = 0, k = 0, j = 0; i < Num_sites; i++) {
    svd_predicate_push_back_C2<K>(t[i], v, k, site_types, j);
  }

  return svd_predicate_ftC2<Predicate,Return_t,FT,
    ITag,Num_sites>(v, site_types);
}

template<template<class Kernel, class MTag> class Predicate,
	 typename Return_t, class K,
	 class Method_tag, unsigned int Num_sites>
Return_t
svd_predicate_C2(const typename K::Site_2 t[])
{
  typedef typename K::FT                 FT;
  typedef typename K::Intersections_tag  ITag;

  FT v[Num_sites * 12];
  char site_types[Num_sites * 2];

  for (unsigned int i = 0, k = 0, j = 0; i < Num_sites; i++) {
    svd_predicate_push_back_C2<K>(t[i], v, k, site_types, j);
  }

  return svd_predicate_ftC2<Predicate,Return_t,FT,
    Method_tag,ITag,Num_sites>(v, site_types);
}

template<template<class Kernel, class MTag> class Predicate,
	 typename Return_t, class K,
	 class Method_tag, typename Data, unsigned int Num_sites>
Return_t
svd_predicate_C2(const typename K::Site_2 t[], Data data)
{
  typedef typename K::FT                 FT;
  typedef typename K::Intersections_tag  ITag;

  FT v[Num_sites * 12];
  char site_types[Num_sites * 2];

  for (unsigned int i = 0, k = 0, j = 0; i < Num_sites; i++) {
    svd_predicate_push_back_C2<K>(t[i], v, k, site_types, j);
  }

  return svd_predicate_ftC2<Predicate,Return_t,FT,
    Method_tag,ITag,Data,Num_sites>(v, site_types, data);
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------


template<class K>
inline
bool svd_are_same_points_C2(const typename K::Site_2& p,
			    const typename K::Site_2& q)
{
  typename K::Site_2 site_vec[] = {p, q};
  return svd_predicate_C2<Svd_are_same_points_C2,bool,K,2>(site_vec);
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

template<class K>
inline
Orientation svd_orientation_C2(const typename K::Site_2& p,
			       const typename K::Site_2& q,
			       const typename K::Site_2& r)
{
  typename K::Site_2 site_vec[] = {p, q, r};
  return svd_predicate_C2<Svd_orientation_C2,Orientation,K,3>(site_vec);
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

template<class K, class Method_tag>
inline
Oriented_side
svd_oriented_side_of_bisector_ftC2(const typename K::Site_2& p,
				   const typename K::Site_2& q,
				   const typename K::Site_2& t,
				   Method_tag mtag)
{
  typename K::Site_2 site_vec[] = {p, q, t};
  return svd_predicate_C2<Svd_oriented_side_of_bisector_C2,
    Oriented_side,K,Method_tag,3>(site_vec);
}


//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

template<class K, class Method_tag>
inline
Sign
svd_vertex_conflict_ftC2(const typename K::Site_2& p,
			 const typename K::Site_2& q,
			 const typename K::Site_2& t,
			 Method_tag mtag)
{
  typename K::Site_2 site_vec[] = {p, q, t};
  return
    svd_predicate_C2<Svd_incircle_2,Sign,K,Method_tag,3>(site_vec);
}

//--------------------------------------------------------------------------

template<class K, class Method_tag>
inline
Sign
svd_vertex_conflict_ftC2(const typename K::Site_2& p,
			 const typename K::Site_2& q,
			 const typename K::Site_2& r,
			 const typename K::Site_2& t,
			 Method_tag mtag)
{
  typename K::Site_2 site_vec[] = {p, q, r, t};
  return
    svd_predicate_C2<Svd_incircle_2,Sign,K,Method_tag,4>(site_vec);
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

template<class K, class Method_tag>
inline
bool
svd_finite_edge_conflict_ftC2(const typename K::Site_2& p,
			      const typename K::Site_2& q,
			      const typename K::Site_2& t,
			      Sign sgn, Method_tag mtag)
{
  typename K::Site_2 site_vec[] = {p, q, t};
  return svd_predicate_C2<Svd_finite_edge_interior_2,bool,K,
    Method_tag,Sign,3>(site_vec, sgn);
}

//--------------------------------------------------------------------------

template<class K, class Method_tag>
inline
bool
svd_finite_edge_conflict_ftC2(const typename K::Site_2& p,
			      const typename K::Site_2& q,
			      const typename K::Site_2& r,
			      const typename K::Site_2& t,
			      Sign sgn, Method_tag mtag)
{
  typename K::Site_2 site_vec[] = {p, q, r, t};
  return svd_predicate_C2<Svd_finite_edge_interior_2,bool,K,
    Method_tag,Sign,4>(site_vec, sgn);
}


//--------------------------------------------------------------------------

template<class K, class Method_tag>
inline
bool
svd_finite_edge_conflict_ftC2(const typename K::Site_2& p,
			      const typename K::Site_2& q,
			      const typename K::Site_2& r,
			      const typename K::Site_2& s,
			      const typename K::Site_2& t,
			      Sign sgn, Method_tag mtag)
{
  typename K::Site_2 site_vec[] = {p, q, r, s, t};
  return svd_predicate_C2<Svd_finite_edge_interior_2,bool,K,
    Method_tag,Sign,5>(site_vec, sgn);
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

template<class K, class Method_tag>
inline
bool
svd_infinite_edge_conflict_ftC2(const typename K::Site_2& q,
				const typename K::Site_2& r,
				const typename K::Site_2& s,
				const typename K::Site_2& t,
				Sign sgn, Method_tag mtag)
{
  typename K::Site_2 site_vec[] = {q, r, s, t};
  return svd_predicate_C2<Svd_infinite_edge_interior_2,bool,K,
    Method_tag,Sign,4>(site_vec, sgn);
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

template<class K, class Method_tag>
inline
bool
svd_is_degenerate_edge_ftC2(const typename K::Site_2& p,
			    const typename K::Site_2& q,
			    const typename K::Site_2& r,
			    const typename K::Site_2& t,
			    Method_tag mtag)
{
  typename K::Site_2 site_vec[] = {p, q, r, t};
  return
    svd_predicate_C2<Svd_is_degenerate_edge_C2,bool,K,Method_tag,4>(site_vec);
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

template<class K, class Method_tag>
inline
std::pair<int,int>
svd_arrangement_type_C2(const typename K::Site_2& p,
			const typename K::Site_2& q,
			Method_tag mtag)
{
  typename K::Site_2 site_vec[2] = {p, q};
  return
    svd_predicate_C2<Svd_arrangement_type_C2,std::pair<int,int>,K,2>(site_vec);
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

template<class K>
inline
bool
svd_are_parallel_C2(const typename K::Site_2& p,
		    const typename K::Site_2& q)
{
#if 1
  typename K::Site_2 site_vec[2] = {p, q};
  return
    svd_predicate_C2<Svd_are_parallel_C2,bool,K,2>(site_vec);

#else
  typedef typename K::Segment_2  Segment_2;

  CGAL_precondition( p.is_segment() && q.is_segment() );

  Segment_2 s1 = p.segment();
  Segment_2 s2 = q.segment();

  return svd_are_parallel_ftC2(s1.source().x(),	s1.source().y(),
			       s1.target().x(),	s1.target().y(),
			       s2.source().x(),	s2.source().y(),
			       s2.target().x(),	s2.target().y());
#endif
}


template<class FT>
inline
bool
svd_are_parallel_ftC2(const FT& x1, const FT& y1,
		      const FT& x2, const FT& y2,
		      const FT& x3, const FT& y3,
		      const FT& x4, const FT& y4)
{
  must_be_filtered(FT());

  FT det = det2x2_by_formula(x2 - x1, x4 - x3,
			     y2 - y1, y4 - y3);

  return ( CGAL::sign(det) == CGAL::ZERO );
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

template<class K, class Method_tag>
inline
Oriented_side
svd_oriented_side_ftC2(const typename K::Site_2& s1,
		       const typename K::Site_2& s2,
		       const typename K::Site_2& s3,
		       const typename K::Site_2& s,
		       const typename K::Site_2& p,
		       Method_tag mtag)
{
  typename K::Site_2 site_vec[] = {s1, s2, s3, s, p};
  return svd_predicate_C2<Svd_oriented_side_C2,Oriented_side,K,
    Method_tag,5>(site_vec);
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------


CGAL_END_NAMESPACE


#ifdef CGAL_ARITHMETIC_FILTER_H
#ifndef CGAL_ARITHMETIC_FILTER_SVD_PREDICATES_FTC2_H
#include <CGAL/Arithmetic_filter/predicates/svd_predicates_ftC2.h>
#endif // CGAL_ARITHMETIC_FILTER_SVD_PREDICATES_FTC2_H
#endif

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_PREDICATES_FTC2_H


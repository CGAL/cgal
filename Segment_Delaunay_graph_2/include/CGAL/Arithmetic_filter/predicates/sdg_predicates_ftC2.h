// Copyright (c) 2003,2004,2005  INRIA Sophia-Antipolis (France) and
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




#ifndef CGAL_ARITHMETIC_FILTER_SDG_PREDICATES_FTC2_H
#define CGAL_ARITHMETIC_FILTER_SDG_PREDICATES_FTC2_H


#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Filtered_exact.h>

#include <CGAL/Segment_Delaunay_graph_2/basic.h>
#include <CGAL/Segment_Delaunay_graph_2/Predicates_ftC2.h>

CGAL_BEGIN_NAMESPACE

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

CGAL_SEGMENT_DELAUNAY_GRAPH_2_BEGIN_NAMESPACE

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES
static unsigned int num_failures_are_same_points = 0;
static unsigned int num_failures_side_of_bisector = 0;
static unsigned int num_failures_vertex_conflict = 0;
static unsigned int num_failures_finite_edge_conflict = 0;
static unsigned int num_failures_infinite_edge_conflict = 0;
static unsigned int num_failures_is_degenerate_edge = 0;
static unsigned int num_failures_arrangement_type = 0;
static unsigned int num_failures_are_parallel = 0;
static unsigned int num_failures_oriented_side = 0;
#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES

CGAL_SEGMENT_DELAUNAY_GRAPH_2_END_NAMESPACE

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------

template < class CT, class ET, bool Protected, class Cache, class ITag>
class Sdg_are_same_points_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Sdg_are_same_points_ftC2<IT,ITag>                 IT_Predicate;
  typedef Sdg_are_same_points_ftC2<ET,ITag>                 ET_Predicate;


public:
  inline
  bool operator()(const FT v[], const char site_types[]) const
  {
    try {
      Protect_FPU_rounding<Protected> Protection;

      IT v_IT[24];

      for (unsigned int i = 0; i < 24; i++) {
	v_IT[i] = v[i].interval();
      }
      return IT_Predicate()(v_IT, site_types);
    }
    catch (Interval_nt_advanced::unsafe_comparison) {
      Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES
      CGAL_SEGMENT_DELAUNAY_GRAPH_2_INTERNAL_NS::num_failures_are_same_points++;
#endif

      ET v_ET[24];

      for (unsigned int i = 0; i < 24; i++) {
	v_ET[i] = v[i].exact();
      }

      return ET_Predicate()(v_ET, site_types);
    }
  }
};

//----------------------------------------------------------------------------

template < class CT, class ET, bool Protected, class Cache, class ITag>
class Sdg_orientation_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Sdg_orientation_ftC2<IT,ITag>                     IT_Predicate;
  typedef Sdg_orientation_ftC2<ET,ITag>                     ET_Predicate;


public:
  inline
  Orientation operator()(const FT v[], const char site_types[]) const
  {
    try {
      Protect_FPU_rounding<Protected> Protection;

      IT v_IT[36];

      for (unsigned int i = 0; i < 36; i++) {
	v_IT[i] = v[i].interval();
      }
      return IT_Predicate()(v_IT, site_types);
    }
    catch (Interval_nt_advanced::unsafe_comparison) {
      Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES
      CGAL_SEGMENT_DELAUNAY_GRAPH_2_INTERNAL_NS::num_failures_are_same_points++;
#endif

      ET v_ET[36];

      for (unsigned int i = 0; i < 36; i++) {
	v_ET[i] = v[i].exact();
      }

      return ET_Predicate()(v_ET, site_types);
    }
  }
};

//----------------------------------------------------------------------------

template < class CT, class ET, bool Protected, class Cache,
	   class MTag, class ITag>
class Sdg_oriented_side_of_bisector_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Sdg_oriented_side_of_bisector_ftC2<IT,Sqrt_field_tag,ITag>
  IT_Predicate;

  typedef Sdg_oriented_side_of_bisector_ftC2<ET,MTag,ITag>  ET_Predicate;

public:
  inline
  Oriented_side operator()(const FT v[], const char site_types[]) const
  {
    try {
      Protect_FPU_rounding<Protected> Protection;

      IT v_IT[36];

      for (unsigned int i = 0; i < 36; i++) {
	v_IT[i] = v[i].interval();
      }
      return IT_Predicate()(v_IT, site_types);
    }
    catch (Interval_nt_advanced::unsafe_comparison) {
      Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES
      CGAL_SEGMENT_DELAUNAY_GRAPH_2_INTERNAL_NS::num_failures_side_of_bisector++;
#endif

      ET v_ET[36];

      for (unsigned int i = 0; i < 36; i++) {
	v_ET[i] = v[i].exact();
      }

      return ET_Predicate()(v_ET, site_types);
    }
  }
};

//----------------------------------------------------------------------------

template < class CT, class ET, bool Protected, class Cache, class MTag,
	   class ITag, int Num_sites>
class Sdg_vertex_conflict_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag,Num_sites>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Sdg_vertex_conflict_ftC2<IT,Sqrt_field_tag,ITag,Num_sites>
  IT_Predicate;

  typedef Sdg_vertex_conflict_ftC2<ET,MTag,ITag,Num_sites>  ET_Predicate;

public:
  inline
  Sign operator()(const FT v[], const char site_types[]) const
  {
    try {
      Protect_FPU_rounding<Protected> Protection;

      IT v_IT[12 * Num_sites];

      for (unsigned int i = 0; i < 12 * Num_sites; i++) {
	v_IT[i] = v[i].interval();
      }

      return IT_Predicate()(v_IT, site_types);
    }
    catch (Interval_nt_advanced::unsafe_comparison) {
      Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES
      CGAL_SEGMENT_DELAUNAY_GRAPH_2_INTERNAL_NS::num_failures_vertex_conflict++;
#endif

      ET v_ET[12 * Num_sites];

      for (unsigned int i = 0; i < 12 * Num_sites; i++) {
	v_ET[i] = v[i].exact();
      }

      return ET_Predicate()(v_ET, site_types);
    }
  }
};

//----------------------------------------------------------------------------


template < class CT, class ET, bool Protected, class Cache,
	   class MTag, class ITag, int Num_sites>
class Sdg_finite_edge_conflict_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag,Num_sites>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Sdg_finite_edge_conflict_ftC2<IT,Sqrt_field_tag,ITag,Num_sites>
  IT_Predicate;

  typedef Sdg_finite_edge_conflict_ftC2<ET,MTag,ITag,Num_sites>  ET_Predicate;

public:
  inline
  bool operator()(const FT v[],	Sign sgn, char site_types[]) const
  {
    try {
      Protect_FPU_rounding<Protected> Protection;

      IT v_IT[12 * Num_sites];

      for (unsigned int i = 0; i < 12 * Num_sites; i++) {
	v_IT[i] = v[i].interval();
      }
      return IT_Predicate()(v_IT, sgn, site_types);
    }
    catch (Interval_nt_advanced::unsafe_comparison) {
      Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES
      CGAL_SEGMENT_DELAUNAY_GRAPH_2_INTERNAL_NS::num_failures_finite_edge_conflict++;
#endif

      ET v_ET[12 * Num_sites];

      for (unsigned int i = 0; i < 12 * Num_sites; i++) {
	v_ET[i] = v[i].exact();
      }

      return ET_Predicate()(v_ET, sgn, site_types);
    }
  }
};


//----------------------------------------------------------------------------


template <class CT, class ET, bool Protected, class Cache, class MTag,
	  class ITag>
class Sdg_infinite_edge_conflict_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Sdg_infinite_edge_conflict_ftC2<IT,Sqrt_field_tag,ITag>
  IT_Predicate;

  typedef Sdg_infinite_edge_conflict_ftC2<ET,MTag,ITag>     ET_Predicate;

public:
  inline
  bool operator()(const FT v[], Sign sgn, const char site_types[]) const
  {
    try {
      Protect_FPU_rounding<Protected> Protection;

      IT v_IT[48];

      for (unsigned int i = 0; i < 48; i++) {
	v_IT[i] = v[i].interval();
      }
      return IT_Predicate()(v_IT, sgn, site_types);
    }
    catch (Interval_nt_advanced::unsafe_comparison) {
      Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES
      CGAL_SEGMENT_DELAUNAY_GRAPH_2_INTERNAL_NS::num_failures_infinite_edge_conflict++;
#endif

      ET v_ET[48];

      for (unsigned int i = 0; i < 48; i++) {
	v_ET[i] = v[i].exact();
      }

      return ET_Predicate()(v_ET, sgn, site_types);
    }
  }
};


//----------------------------------------------------------------------------

template <class CT, class ET, bool Protected, class Cache, class MTag,
	  class ITag>
class Sdg_is_degenerate_edge_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Sdg_is_degenerate_edge_ftC2<IT,Sqrt_field_tag,ITag>  IT_Predicate;
  typedef Sdg_is_degenerate_edge_ftC2<ET,MTag,ITag>            ET_Predicate;

public:
  inline
  bool operator()(const FT v[], const char site_types[]) const 
  {
    try {
      Protect_FPU_rounding<Protected> Protection;

      IT v_IT[48];

      for (unsigned int i = 0; i < 48; i++) {
	v_IT[i] = v[i].interval();
      }
      return IT_Predicate()(v_IT, site_types);
    }
    catch (Interval_nt_advanced::unsafe_comparison) {
      Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES
      CGAL_SEGMENT_DELAUNAY_GRAPH_2_INTERNAL_NS::num_failures_is_degenerate_edge++;
#endif

      ET v_ET[48];

      for (unsigned int i = 0; i < 48; i++) {
	v_ET[i] = v[i].exact();
      }

      return ET_Predicate()(v_ET, site_types);
    }
  }
};

//----------------------------------------------------------------------------

template <class CT, class ET, bool Protected, class Cache, class MTag,
	  class ITag>
class Sdg_arrangement_type_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Sdg_arrangement_type_ftC2<IT,Sqrt_field_tag,ITag>  IT_Predicate;
  typedef Sdg_arrangement_type_ftC2<ET,MTag,ITag>            ET_Predicate;

public:
  typedef typename IT_Predicate::result_type   result_type;

  inline
  result_type operator()(const FT v[],  const char site_types[]) const
  {
    try {
      Protect_FPU_rounding<Protected> Protection;

      IT v_IT[24];

      for (unsigned int i = 0; i < 24; i++) {
	v_IT[i] = v[i].interval();
      }
      return IT_Predicate()(v_IT, site_types);
    }
    catch (Interval_nt_advanced::unsafe_comparison) {
      Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES
      CGAL_SEGMENT_DELAUNAY_GRAPH_2_INTERNAL_NS::num_failures_arrangement_type++;
#endif

      ET v_ET[24];

      for (unsigned int i = 0; i < 24; i++) {
	v_ET[i] = v[i].exact();
      }

      return ET_Predicate()(v_ET, site_types);
    }
  }
};

//----------------------------------------------------------------------------

template <class CT, class ET, bool Protected, class Cache, class ITag>
class Sdg_are_parallel_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Sdg_are_parallel_ftC2<IT,ITag>                    IT_Predicate;
  typedef Sdg_are_parallel_ftC2<ET,ITag>                    ET_Predicate;

public:
  inline
  bool operator()(const FT v[],	const char site_types[]) const
  {
    try {
      Protect_FPU_rounding<Protected> Protection;

      IT v_IT[24];

      for (unsigned int i = 0; i < 24; i++) {
	v_IT[i] = v[i].interval();
      }
      return IT_Predicate()(v_IT, site_types);
    }
    catch (Interval_nt_advanced::unsafe_comparison) {
      Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES
      CGAL_SEGMENT_DELAUNAY_GRAPH_2_INTERNAL_NS::num_failures_are_parallel++;
#endif

      ET v_ET[24];

      for (unsigned int i = 0; i < 24; i++) {
	v_ET[i] = v[i].exact();
      }

      return ET_Predicate()(v_ET, site_types);
    }
  }
};

//----------------------------------------------------------------------------

template <class CT, class ET, bool Protected, class Cache, class MTag,
	  class ITag, int Num_sites>
class Sdg_oriented_side_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag,Num_sites>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Sdg_oriented_side_ftC2<IT,Sqrt_field_tag,ITag,Num_sites>
  IT_Predicate;

  typedef Sdg_oriented_side_ftC2<ET,MTag,ITag,Num_sites>    ET_Predicate;

public:
  inline
  Oriented_side operator()(const FT v[], const char site_types[]) const
  {
    try {
      Protect_FPU_rounding<Protected> Protection;

      IT v_IT[12 * Num_sites];

      for (unsigned int i = 0; i < 12 * Num_sites; i++) {
	v_IT[i] = v[i].interval();
      }
      return IT_Predicate()(v_IT, site_types);
    }
    catch (Interval_nt_advanced::unsafe_comparison) {
      Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTER_FAILURES
      CGAL_SEGMENT_DELAUNAY_GRAPH_2_INTERNAL_NS::num_failures_oriented_side++;
#endif

      ET v_ET[12 * Num_sites];

      for (unsigned int i = 0; i < 12 * Num_sites; i++) {
	v_ET[i] = v[i].exact();
      }

      return ET_Predicate()(v_ET, site_types);
    }
  }
};

//----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_ARITHMETIC_FILTER_SDG_PREDICATES_FTC2_H

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




#ifndef CGAL_ARITHMETIC_FILTER_SVD_PREDICATES_FTC2_H
#define CGAL_ARITHMETIC_FILTER_SVD_PREDICATES_FTC2_H


#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Filtered_exact.h>

#include <CGAL/predicates/Segment_Voronoi_diagram_predicates_ftC2.h>

CGAL_BEGIN_NAMESPACE

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

static unsigned int num_failures_are_same_points = 0;
static unsigned int num_failures_side_of_bisector = 0;
static unsigned int num_failures_vertex_conflict = 0;
static unsigned int num_failures_finite_edge_conflict = 0;
static unsigned int num_failures_infinite_edge_conflict = 0;
static unsigned int num_failures_is_degenerate_edge = 0;
static unsigned int num_failures_arrangement_type = 0;
static unsigned int num_failures_are_parallel = 0;
static unsigned int num_failures_oriented_side = 0;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------

template < class CT, class ET, bool Protected, class Cache, class ITag>
class Svd_are_same_points_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Svd_are_same_points_ftC2<IT,ITag>                 IT_Predicate;
  typedef Svd_are_same_points_ftC2<ET,ITag>                 ET_Predicate;


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

      num_failures_are_same_points++;

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
class Svd_orientation_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Svd_orientation_ftC2<IT,ITag>                     IT_Predicate;
  typedef Svd_orientation_ftC2<ET,ITag>                     ET_Predicate;


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

      num_failures_are_same_points++;

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
class Svd_oriented_side_of_bisector_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Svd_oriented_side_of_bisector_ftC2<IT,Sqrt_field_tag,ITag>
  IT_Predicate;

  typedef Svd_oriented_side_of_bisector_ftC2<ET,MTag,ITag>  ET_Predicate;

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

      num_failures_side_of_bisector++;

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
class Svd_vertex_conflict_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag,Num_sites>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Svd_vertex_conflict_ftC2<IT,Sqrt_field_tag,ITag,Num_sites>
  IT_Predicate;

  typedef Svd_vertex_conflict_ftC2<ET,MTag,ITag,Num_sites>  ET_Predicate;

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

      num_failures_vertex_conflict++;

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
class Svd_finite_edge_conflict_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag,Num_sites>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Svd_finite_edge_conflict_ftC2<IT,Sqrt_field_tag,ITag,Num_sites>
  IT_Predicate;

  typedef Svd_finite_edge_conflict_ftC2<ET,MTag,ITag,Num_sites>  ET_Predicate;

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

      num_failures_finite_edge_conflict++;

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
class Svd_infinite_edge_conflict_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Svd_infinite_edge_conflict_ftC2<IT,Sqrt_field_tag,ITag>
  IT_Predicate;

  typedef Svd_infinite_edge_conflict_ftC2<ET,MTag,ITag>     ET_Predicate;

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

      num_failures_infinite_edge_conflict++;

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
class Svd_is_degenerate_edge_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Svd_is_degenerate_edge_ftC2<IT,Sqrt_field_tag,ITag>  IT_Predicate;
  typedef Svd_is_degenerate_edge_ftC2<ET,MTag,ITag>            ET_Predicate;

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

      num_failures_is_degenerate_edge++;

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
class Svd_arrangement_type_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Svd_arrangement_type_ftC2<IT,Sqrt_field_tag,ITag>  IT_Predicate;
  typedef Svd_arrangement_type_ftC2<ET,MTag,ITag>            ET_Predicate;

public:

  inline
  std::pair<int,int>
  operator()(const FT v[],  const char site_types[]) const
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

      num_failures_arrangement_type++;

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
class Svd_are_parallel_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Svd_are_parallel_ftC2<IT,ITag>                    IT_Predicate;
  typedef Svd_are_parallel_ftC2<ET,ITag>                    ET_Predicate;

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

      num_failures_are_parallel++;

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
	  class ITag>
class Svd_oriented_side_ftC2<
  Filtered_exact<CT,ET,Dynamic,Protected,Cache>,
  MTag,ITag>
{
private:
  typedef Interval_nt_advanced                              IT;
  typedef Filtered_exact<CT,ET,Dynamic,Protected,Cache>     FT;

  typedef Svd_oriented_side_ftC2<IT,Sqrt_field_tag,ITag>    IT_Predicate;
  typedef Svd_oriented_side_ftC2<ET,MTag,ITag>              ET_Predicate;

public:
  inline
  Oriented_side operator()(const FT v[], const char site_types[]) const
  {
    try {
      Protect_FPU_rounding<Protected> Protection;

      IT v_IT[60];

      for (unsigned int i = 0; i < 60; i++) {
	v_IT[i] = v[i].interval();
      }
      return IT_Predicate()(v_IT, site_types, mtag, itag);
    }
    catch (Interval_nt_advanced::unsafe_comparison) {
      Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

      num_failures_oriented_side++;

      ET v_ET[60];

      for (unsigned int i = 0; i < 60; i++) {
	v_ET[i] = v[i].exact();
      }

      return ET_Predicate()(v_ET, site_types, mtag, itag);
    }
  }
};

//----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_ARITHMETIC_FILTER_SVD_PREDICATES_FTC2_H

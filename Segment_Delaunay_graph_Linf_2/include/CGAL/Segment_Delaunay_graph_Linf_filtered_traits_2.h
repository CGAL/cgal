// Copyright (c) 2015  Università della Svizzera italiana.
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_FILTERED_TRAITS_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_FILTERED_TRAITS_2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Segment_Delaunay_graph_Linf_2/basic.h>

#include <CGAL/Segment_Delaunay_graph_Linf_2/Filtered_traits_base_2.h>
#include <CGAL/Segment_Delaunay_graph_2/Filtered_traits_concept_check_tags.h>

// includes for the default parameters of the filtered traits
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>
#else
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/number_utils_classes.h>

namespace CGAL {

#define SDG2_INS CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Internal

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// the filtered Traits classes
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

// this traits class does support intersecting segments
template<class CK,
	 class CK_MTag = Field_with_sqrt_tag,
#ifdef CGAL_USE_GMP
	 class EK      = Simple_cartesian< Gmpq >,
#else
	 class EK      = Simple_cartesian< Quotient<MP_Float> >,
#endif
	 class EK_MTag = Field_tag,
	 class FK      = Simple_cartesian< Interval_nt<false> >,
	 class FK_MTag = Field_with_sqrt_tag,
	 class C2E     = Cartesian_converter<CK, EK>,
	 class C2F     =
	 Cartesian_converter<CK, FK, To_interval<typename CK::RT> > >
struct Segment_Delaunay_graph_Linf_filtered_traits_2
  : public Segment_Delaunay_graph_Linf_filtered_traits_base_2<CK, CK_MTag,
							 EK, EK_MTag,
							 FK, FK_MTag,
							 C2E, C2F,
							 Tag_true>
{
public:
  Segment_Delaunay_graph_Linf_filtered_traits_2() {
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,CK_MTag,2>();
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,EK_MTag,4>();
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,FK_MTag,6>();
  }
};


template<class CK, class EK, class EK_MTag, class FK, class FK_MTag,
	 class C2E, class C2F>
struct Segment_Delaunay_graph_Linf_filtered_traits_2<CK, Field_tag,
						EK, EK_MTag,
						FK, FK_MTag,
						C2E, C2F>
  : public Segment_Delaunay_graph_Linf_filtered_traits_base_2<CK, Integral_domain_without_division_tag,
							 EK, EK_MTag,
							 FK, FK_MTag,
							 C2E, C2F,
							 Tag_true>
{
public:
  Segment_Delaunay_graph_Linf_filtered_traits_2() {
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,EK_MTag,4>();
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,FK_MTag,6>();
  }
};

template<class CK, class CK_MTag, class EK, class FK, class FK_MTag,
	 class C2E, class C2F>
struct Segment_Delaunay_graph_Linf_filtered_traits_2<CK, CK_MTag,
						EK, Field_tag,
						FK, FK_MTag,
						C2E, C2F>
  : public Segment_Delaunay_graph_Linf_filtered_traits_base_2<CK, CK_MTag,
							 EK, Integral_domain_without_division_tag,
							 FK, FK_MTag,
							 C2E, C2F,
							 Tag_true>
{
public:
  Segment_Delaunay_graph_Linf_filtered_traits_2() {
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,CK_MTag,2>();
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,FK_MTag,6>();
  }
};

template<class CK, class CK_MTag, class EK, class EK_MTag, class FK,
	 class C2E, class C2F>
struct Segment_Delaunay_graph_Linf_filtered_traits_2<CK, CK_MTag,
						EK, EK_MTag,
						FK, Field_tag,
						C2E, C2F>
  : public Segment_Delaunay_graph_Linf_filtered_traits_base_2<CK, CK_MTag,
							 EK, EK_MTag,
							 FK, Integral_domain_without_division_tag,
							 C2E, C2F,
							 Tag_true>
{
public:
  Segment_Delaunay_graph_Linf_filtered_traits_2() {
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,CK_MTag,2>();
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,EK_MTag,4>();
  }
};

template<class CK, class CK_MTag, class EK, class FK,
	 class C2E, class C2F>
struct Segment_Delaunay_graph_Linf_filtered_traits_2<CK, CK_MTag,
						EK, Field_tag,
						FK, Field_tag,
						C2E, C2F>
  : public Segment_Delaunay_graph_Linf_filtered_traits_base_2<CK, CK_MTag,
							 EK, Integral_domain_without_division_tag,
							 FK, Integral_domain_without_division_tag,
							 C2E, C2F,
							 Tag_true>
{
public:
  Segment_Delaunay_graph_Linf_filtered_traits_2() {
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,CK_MTag,2>();
  }
};

template<class CK, class EK, class EK_MTag, class FK,
	 class C2E, class C2F>
struct Segment_Delaunay_graph_Linf_filtered_traits_2<CK, Field_tag,
						 EK, EK_MTag,
						 FK, Field_tag,
						 C2E, C2F>
  : public Segment_Delaunay_graph_Linf_filtered_traits_base_2<CK, Integral_domain_without_division_tag,
							 EK, EK_MTag,
							 FK, Integral_domain_without_division_tag,
							 C2E, C2F,
							 Tag_true>
{
public:
  Segment_Delaunay_graph_Linf_filtered_traits_2() {
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,EK_MTag,4>();
  }
};

template<class CK, class EK, class FK, class FK_MTag,
	 class C2E, class C2F>
struct Segment_Delaunay_graph_Linf_filtered_traits_2<CK, Field_tag,
						EK, Field_tag,
						FK, FK_MTag,
						C2E, C2F>
  : public Segment_Delaunay_graph_Linf_filtered_traits_base_2<CK, Integral_domain_without_division_tag,
							 EK, Integral_domain_without_division_tag,
							 FK, FK_MTag,
							 C2E, C2F,
							 Tag_true>
{
public:
  Segment_Delaunay_graph_Linf_filtered_traits_2() {
    SDG2_INS::Concept_check_tags<Integral_domain_without_division_tag,FK_MTag,6>();
  }
};

template<class CK, class EK, class FK, class C2E, class C2F>
struct Segment_Delaunay_graph_Linf_filtered_traits_2<CK, Field_tag,
						EK, Field_tag,
						FK, Field_tag,
						C2E, C2F>
  : public Segment_Delaunay_graph_Linf_filtered_traits_base_2<CK, Integral_domain_without_division_tag,
							 EK, Integral_domain_without_division_tag,
							 FK, Integral_domain_without_division_tag,
							 C2E, C2F,
							 Tag_true>
{};

//=========================================================================


// this traits class does NOT support intersecting segments
template<class CK,
	 class CK_MTag = Field_with_sqrt_tag,
#ifdef CGAL_USE_GMP
	 class EK      = Simple_cartesian< Gmpq >,
#else
	 class EK      = Simple_cartesian< MP_Float >,
#endif
	 class EK_MTag = Integral_domain_without_division_tag,
	 class FK      = Simple_cartesian< Interval_nt<false> >,
	 class FK_MTag = Field_with_sqrt_tag,
	 class C2E     = Cartesian_converter<CK, EK>,
	 class C2F     =
	 Cartesian_converter<CK, FK, To_interval<typename CK::RT> > >
struct Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2
  : public Segment_Delaunay_graph_Linf_filtered_traits_base_2<CK, CK_MTag,
							 EK, EK_MTag,
							 FK, FK_MTag,
							 C2E, C2F,
							 Tag_false>
{
  Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2() {
    SDG2_INS::Concept_check_tags_wi<Integral_domain_without_division_tag,CK_MTag,2>();
    SDG2_INS::Concept_check_tags_wi<Integral_domain_without_division_tag,EK_MTag,4>();
    SDG2_INS::Concept_check_tags_wi<Integral_domain_without_division_tag,FK_MTag,6>();
  }
};

#undef SDG2_INS

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_FILTERED_TRAITS_2_H

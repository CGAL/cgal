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



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_FILTERED_TRAITS_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_FILTERED_TRAITS_2_H

#include <CGAL/Segment_Voronoi_diagram_short_names_2.h>

#include <CGAL/Segment_Voronoi_diagram_filtered_traits_base_2.h>

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


CGAL_BEGIN_NAMESPACE

#if defined(__sun) && defined(__SUNPRO_CC)
// workaround for the Sun CC-5.30 compiler; it does not like default
// template parameters that are themselves templates and have
// templated classes as parameters, which have then nested types as
// arguments... oooof!!!
//
// In case you did understand what I just described you are most
// probably crazy... If you did not, look below to see what kind of
// code CC-5.30 did not like.
namespace CGALi {

  template<class CK, class FK>
  struct SVD_SUNPRO_CC_Interval_converter
    : public Cartesian_converter<CK, FK,
                                 To_interval< typename CK::RT > >
  {
  };

}
#endif


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
	 class CK_MTag = Sqrt_field_tag,
#ifdef CGAL_USE_GMP
	 class EK      = Simple_cartesian< Gmpq >,
#else
	 class EK      = Simple_cartesian< Quotient<MP_Float> >,
#endif
	 class EK_MTag = Ring_tag,
	 class FK      = Simple_cartesian< Interval_nt<false> >,
	 class FK_MTag = Sqrt_field_tag,
	 class C2E     = Cartesian_converter<CK, EK>,
#if defined(__sun) && defined(__SUNPRO_CC)
         class C2F     = CGALi::SVD_SUNPRO_CC_Interval_converter<CK, FK> > >
#else
	 class C2F     =
	 Cartesian_converter<CK, FK, To_interval<typename CK::RT> > >
#endif
class Segment_Voronoi_diagram_filtered_traits_2
  : public Segment_Voronoi_diagram_filtered_traits_base_2<CK, CK_MTag,
							  EK, EK_MTag,
							  FK, FK_MTag,
							  C2E, C2F,
							  Tag_true>
{};


// this traits class does NOT support intersecting segments
template<class CK,
	 class CK_MTag = Sqrt_field_tag,
#ifdef CGAL_USE_GMP
	 class EK      = Simple_cartesian< Gmpq >,
#else
	 class EK      = Simple_cartesian< MP_Float >,
#endif
	 class EK_MTag = Ring_tag,
	 class FK      = Simple_cartesian< Interval_nt<false> >,
	 class FK_MTag = Sqrt_field_tag,
	 class C2E     = Cartesian_converter<CK, EK>,
#if defined(__sun) && defined(__SUNPRO_CC)
         class C2F     = CGALi::SVD_SUNPRO_CC_Interval_converter<CK, FK> > >
#else
	 class C2F     =
	 Cartesian_converter<CK, FK, To_interval<typename CK::RT> > >
#endif
class Segment_Voronoi_diagram_filtered_traits_without_intersections_2
  : public Segment_Voronoi_diagram_filtered_traits_base_2<CK, CK_MTag,
							  EK, EK_MTag,
							  FK, FK_MTag,
							  C2E, C2F,
							  Tag_false>
{};



CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_FILTERED_TRAITS_2_H

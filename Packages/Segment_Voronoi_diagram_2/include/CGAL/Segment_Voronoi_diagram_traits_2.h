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



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_2_H


#include <CGAL/Segment_Voronoi_diagram_short_names_2.h>

#include <CGAL/Number_type_traits.h>
#include <CGAL/Segment_Voronoi_diagram_traits_base_2.h>

#include <CGAL/Filtered_kernel.h>
#include <CGAL/Segment_Voronoi_diagram_filtered_traits_2.h>


CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------
// the Traits classes
//-----------------------------------------------------------------------

// this traits class does support intersecting segments
template<class R, class MTag = Field_tag>
struct Segment_Voronoi_diagram_traits_2
  : public Segment_Voronoi_diagram_traits_base_2<R,MTag,Tag_true> {};

template<class R>
struct Segment_Voronoi_diagram_traits_2<R,Field_tag>
  : public Segment_Voronoi_diagram_traits_base_2<R,Ring_tag,Tag_true> {};


// Concept checking
template<class R>
struct Segment_Voronoi_diagram_traits_2<R,Ring_tag>
  : public Segment_Voronoi_diagram_traits_base_2<R,Ring_tag,Tag_true>
{
  Segment_Voronoi_diagram_traits_2() {
    THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Field_tag_OR_Sqrt_field_tag
    ( R() );
  }
};

template<class R>
struct Segment_Voronoi_diagram_traits_2<R,Euclidean_ring_tag>
  : public Segment_Voronoi_diagram_traits_base_2<R,Ring_tag,Tag_true>
{
  Segment_Voronoi_diagram_traits_2() {
    THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Field_tag_OR_Sqrt_field_tag
    ( R() );
  }
};

// Specializations for filtered_kernel
template<class R>
struct Segment_Voronoi_diagram_traits_2<Filtered_kernel<R>,Field_tag>
  : public
  Segment_Voronoi_diagram_filtered_traits_2<R,Field_tag,
					    typename Filtered_kernel<R>::EK,
					    Field_tag,
					    typename Filtered_kernel<R>::FK,
					    Sqrt_field_tag>
{};

template<class R>
struct Segment_Voronoi_diagram_traits_2<Filtered_kernel<R>,Sqrt_field_tag>
  : public
  Segment_Voronoi_diagram_filtered_traits_2<R,Sqrt_field_tag,
					    typename Filtered_kernel<R>::EK,
					    Field_tag,
					    typename Filtered_kernel<R>::FK,
					    Sqrt_field_tag>
{};

//=========================================================================

// this traits class does NOT support intersecting segments
template<class R, class MTag = Ring_tag>
struct Segment_Voronoi_diagram_traits_without_intersections_2
  : public Segment_Voronoi_diagram_traits_base_2<R,MTag,Tag_false> {};

template<class R>
struct Segment_Voronoi_diagram_traits_without_intersections_2<R,Field_tag>
  : public Segment_Voronoi_diagram_traits_base_2<R,Ring_tag,Tag_false>
{
  Segment_Voronoi_diagram_traits_without_intersections_2() {
    THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Ring_tag_OR_Sqrt_field_tag
    ( R() );
  }
};

template<class R>
struct
Segment_Voronoi_diagram_traits_without_intersections_2<R,Euclidean_ring_tag>
  : public Segment_Voronoi_diagram_traits_base_2<R,Ring_tag,Tag_false>
{
  Segment_Voronoi_diagram_traits_without_intersections_2() {
    THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Ring_tag_OR_Sqrt_field_tag
    ( R() );
  }
};


// Specialization for filtered_kernel
template<class R>
struct
Segment_Voronoi_diagram_traits_without_intersections_2<Filtered_kernel<R>,
						       Ring_tag>
  : public
  Segment_Voronoi_diagram_filtered_traits_without_intersections_2<R,Ring_tag,
					    typename Filtered_kernel<R>::EK,
					    Ring_tag,
					    typename Filtered_kernel<R>::FK,
					    Sqrt_field_tag>
{};

template<class R>
struct
Segment_Voronoi_diagram_traits_without_intersections_2<Filtered_kernel<R>,
						       Sqrt_field_tag>
  : public
  Segment_Voronoi_diagram_filtered_traits_without_intersections_2<R,
					    Sqrt_field_tag,
					    typename Filtered_kernel<R>::EK,
					    Ring_tag,
					    typename Filtered_kernel<R>::FK,
					    Sqrt_field_tag>
{};

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_2_H

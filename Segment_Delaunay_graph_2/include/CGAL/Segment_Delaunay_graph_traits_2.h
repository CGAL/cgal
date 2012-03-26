// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
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
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_TRAITS_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_TRAITS_2_H


#include <CGAL/Segment_Delaunay_graph_2/basic.h>

#include <CGAL/Segment_Delaunay_graph_2/Traits_base_2.h>

#include <CGAL/Filtered_kernel.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>


namespace CGAL {

//-----------------------------------------------------------------------
// the Traits classes
//-----------------------------------------------------------------------

// this traits class does support intersecting segments

// the following definition is okay when giving Field_with_sqrt_tag
// for all other tags I need a specialization
template<class R,
	 class MTag =
	 typename Algebraic_structure_traits<typename R::FT>::Algebraic_category
	 /*Field_tag*/>
struct Segment_Delaunay_graph_traits_2
  : public Segment_Delaunay_graph_traits_base_2<R,MTag,Tag_true> {};

template<class R>
struct Segment_Delaunay_graph_traits_2<R,Field_with_kth_root_tag>
  : public Segment_Delaunay_graph_traits_base_2<R,Field_with_sqrt_tag,Tag_true> {};

template<class R>
struct Segment_Delaunay_graph_traits_2<R,Field_with_root_of_tag>
  : public Segment_Delaunay_graph_traits_base_2<R,Field_with_sqrt_tag,Tag_true> {};

template<class R>
struct Segment_Delaunay_graph_traits_2<R,Field_tag>
  : public Segment_Delaunay_graph_traits_base_2<R,Integral_domain_without_division_tag,Tag_true> {};


// Concept checking
template<class R>
struct Segment_Delaunay_graph_traits_2<R,Integral_domain_without_division_tag>
  : public Segment_Delaunay_graph_traits_base_2<R,Integral_domain_without_division_tag,Tag_true>
{
  Segment_Delaunay_graph_traits_2() {
    THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Field_tag_OR_Field_with_sqrt_tag
    ( R() );
  }
};

template<class R>
struct Segment_Delaunay_graph_traits_2<R,Euclidean_ring_tag>
  : public Segment_Delaunay_graph_traits_base_2<R,Integral_domain_without_division_tag,Tag_true>
{
  Segment_Delaunay_graph_traits_2() {
    THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Field_tag_OR_Field_with_sqrt_tag
    ( R() );
  }
};

// Specializations for filtered_kernel
template<class R>
struct Segment_Delaunay_graph_traits_2<Filtered_kernel<R>,Field_tag>
  : public
  Segment_Delaunay_graph_filtered_traits_2<R,Field_tag,
					   typename Filtered_kernel<R>::Exact_kernel,
					   Field_tag,
					   typename Filtered_kernel<R>::Approximate_kernel,
					   Field_with_sqrt_tag>
{};

template<class R>
struct Segment_Delaunay_graph_traits_2<Filtered_kernel<R>,Field_with_sqrt_tag>
  : public
  Segment_Delaunay_graph_filtered_traits_2<R,Field_with_sqrt_tag,
					   typename Filtered_kernel<R>::Exact_kernel,
					   Field_tag,
					   typename Filtered_kernel<R>::Approximate_kernel,
					   Field_with_sqrt_tag>
{};

//=========================================================================

// this traits class does NOT support intersecting segments
template<class R,
	 class MTag = 
	 typename Algebraic_structure_traits<typename R::FT>::Algebraic_category
	 /*Integral_domain_without_division_tag*/>
struct Segment_Delaunay_graph_traits_without_intersections_2
  : public Segment_Delaunay_graph_traits_base_2<R,MTag,Tag_false> {};

template<class R>
struct Segment_Delaunay_graph_traits_without_intersections_2<R,Field_tag>
  : public Segment_Delaunay_graph_traits_base_2<R,Integral_domain_without_division_tag,Tag_false>
{
  Segment_Delaunay_graph_traits_without_intersections_2() {
    THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Integral_domain_without_division_tag_OR_Field_with_sqrt_tag
    ( R() );
  }
};

template<class R>
struct
Segment_Delaunay_graph_traits_without_intersections_2<R,Euclidean_ring_tag>
  : public Segment_Delaunay_graph_traits_base_2<R,Integral_domain_without_division_tag,Tag_false>
{
  Segment_Delaunay_graph_traits_without_intersections_2() {
    THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Integral_domain_without_division_tag_OR_Field_with_sqrt_tag
    ( R() );
  }
};


// Specialization for filtered_kernel
template<class R>
struct
Segment_Delaunay_graph_traits_without_intersections_2<Filtered_kernel<R>,
						      Integral_domain_without_division_tag>
  : public
  Segment_Delaunay_graph_filtered_traits_without_intersections_2<R,Integral_domain_without_division_tag,
					    typename Filtered_kernel<R>::Exact_kernel,
					    Integral_domain_without_division_tag,
					    typename Filtered_kernel<R>::Approximate_kernel,
					    Field_with_sqrt_tag>
{};

template<class R>
struct
Segment_Delaunay_graph_traits_without_intersections_2<Filtered_kernel<R>,
						      Field_with_sqrt_tag>
  : public
  Segment_Delaunay_graph_filtered_traits_without_intersections_2<R,
					    Field_with_sqrt_tag,
					    typename Filtered_kernel<R>::Exact_kernel,
					    Integral_domain_without_division_tag,
					    typename Filtered_kernel<R>::Approximate_kernel,
					    Field_with_sqrt_tag>
{};

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_TRAITS_2_H

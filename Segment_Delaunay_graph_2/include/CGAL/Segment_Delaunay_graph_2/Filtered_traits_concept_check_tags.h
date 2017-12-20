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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTERED_TRAITS_CONCEPT_CHECK_TAGS_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTERED_TRAITS_CONCEPT_CHECK_TAGS_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>

namespace CGAL {

namespace SegmentDelaunayGraph_2 {

namespace Internal {

  template<class D, class T, int>
  struct Concept_check_tags {};

  template<class D>
  struct Concept_check_tags<D,Integral_domain_without_division_tag,2>
  {
    Concept_check_tags() {
      THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Field_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

  template<class D>
  struct Concept_check_tags<D,Integral_domain_without_division_tag,4> {
    Concept_check_tags() {
      THE_4TH_TEMPLATE_PARAMETER_MUST_EITHER_BE_Field_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

  template<class D>
  struct Concept_check_tags<D,Integral_domain_without_division_tag,6> {
    Concept_check_tags() {
      THE_6TH_TEMPLATE_PARAMETER_MUST_EITHER_BE_Field_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

  //-------------------------------------------------------------------------

  template<class D>
  struct Concept_check_tags<D,Euclidean_ring_tag,2>
  {
    Concept_check_tags() {
      THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Field_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

  template<class D>
  struct Concept_check_tags<D,Euclidean_ring_tag,4> {
    Concept_check_tags() {
      THE_4TH_TEMPLATE_PARAMETER_MUST_EITHER_BE_Field_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

  template<class D>
  struct Concept_check_tags<D,Euclidean_ring_tag,6> {
    Concept_check_tags() {
      THE_6TH_TEMPLATE_PARAMETER_MUST_EITHER_BE_Field_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

  //=========================================================================

  template<class D, class T, int>
  struct Concept_check_tags_wi {};

  template<class D>
  struct Concept_check_tags_wi<D,Field_tag,2>
  {
    Concept_check_tags_wi() {
      THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Integral_domain_without_division_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

  template<class D>
  struct Concept_check_tags_wi<D,Field_tag,4> {
    Concept_check_tags_wi() {
      THE_4TH_TEMPLATE_PARAMETER_MUST_EITHER_BE_Integral_domain_without_division_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

  template<class D>
  struct Concept_check_tags_wi<D,Field_tag,6> {
    Concept_check_tags_wi() {
      THE_6TH_TEMPLATE_PARAMETER_MUST_EITHER_BE_Integral_domain_without_division_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

  //-------------------------------------------------------------------------

  template<class D>
  struct Concept_check_tags_wi<D,Euclidean_ring_tag,2>
  {
    Concept_check_tags_wi() {
      THE_2ND_TEMPLATE_PARAMETER_MUST_EITHER_BE_Integral_domain_without_division_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

  template<class D>
  struct Concept_check_tags_wi<D,Euclidean_ring_tag,4> {
    Concept_check_tags_wi() {
      THE_4TH_TEMPLATE_PARAMETER_MUST_EITHER_BE_Integral_domain_without_division_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

  template<class D>
  struct Concept_check_tags_wi<D,Euclidean_ring_tag,6> {
    Concept_check_tags_wi() {
      THE_6TH_TEMPLATE_PARAMETER_MUST_EITHER_BE_Integral_domain_without_division_tag_OR_Field_with_sqrt_tag
      ( D() );
    }
  };

} // namespace Internal

} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTERED_TRAITS_CONCEPT_CHECK_TAGS_H

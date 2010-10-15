// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France) and
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



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTERED_TRAITS_CONCEPT_CHECK_TAGS_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_FILTERED_TRAITS_CONCEPT_CHECK_TAGS_H

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

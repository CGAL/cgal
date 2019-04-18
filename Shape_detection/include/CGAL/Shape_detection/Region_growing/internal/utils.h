// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_UTILS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_UTILS_H

#include <CGAL/license/Shape_detection.h>

// Boost headers.
#include <boost/mpl/has_xxx.hpp>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

  template<typename GeomTraits> 
  class Default_sqrt {
    
  private:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

  public:
    FT operator()(const FT value) const { 
      
      CGAL_precondition(value >= FT(0));
      return static_cast<FT>(CGAL::sqrt(CGAL::to_double(value)));
    }
  };

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

  // Case: do_not_use_default = false.
  template<typename GeomTraits, 
  bool do_not_use_default = Has_nested_type_Sqrt<GeomTraits>::value>
  class Get_sqrt {
        
  public:
    using Traits = GeomTraits;
    using Sqrt = Default_sqrt<Traits>;

    static Sqrt sqrt_object(const Traits& ) { 
      return Sqrt();
    }
  };

  // Case: do_not_use_default = true.
  template<typename GeomTraits>
  class Get_sqrt<GeomTraits, true> {
        
  public:
    using Traits = GeomTraits;
    using Sqrt = typename Traits::Sqrt;

    static Sqrt sqrt_object(const Traits& traits) { 
      return traits.sqrt_object();
    }
  };

  template<typename FT>
  struct Compare_scores {
    const std::vector<FT>& m_scores;
      
    Compare_scores(const std::vector<FT>& scores) : 
    m_scores(scores) 
    { }

    bool operator()(const std::size_t i, const std::size_t j) const {
        
      CGAL_precondition(i >= 0 && i < m_scores.size());
      CGAL_precondition(j >= 0 && j < m_scores.size());

      return m_scores[i] > m_scores[j];
    }
  };

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_UTILS_H

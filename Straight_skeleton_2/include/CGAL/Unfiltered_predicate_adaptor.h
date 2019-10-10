// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal.
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
// Author(s)     : Sylvain Pion, Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//

#ifndef CGAL_UNFILTERED_PREDICATE_ADAPTOR_H
#define CGAL_UNFILTERED_PREDICATE_ADAPTOR_H

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/basic.h>

namespace CGAL {

template <class CAP>
class Unfiltered_predicate_adaptor
{
  CAP  Certified_approx_predicate;

public:

  typedef typename CAP::result_type  result_type;

  Unfiltered_predicate_adaptor()
  {}

  // These constructors are used for constructive predicates.
  // You should try to avoid constructive predicates, as they will construct
  // the exact values systematically (in the ctor), rather than lazily.
  template <class O>
  Unfiltered_predicate_adaptor(const O &o1)
    : Certified_approx_predicate(o1)
  {}

  template <class O>
  Unfiltered_predicate_adaptor(const O &o1, const O &o2)
    : Certified_approx_predicate(o1, o2)
  {}

  template <class ... A>
  result_type
  operator()(const A&& ... a) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    return static_cast<result_type>(Certified_approx_predicate(std::forward(a)...));
  }
#endif

  // Idem for more than 9 arguments.  Do it on demand.
};

#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
template <class CAP>
  template <class ... A>
typename Unfiltered_predicate_adaptor<CAP>::result_type
Unfiltered_predicate_adaptor<CAP>::
  operator()(const A&& ... a) const
{
  return static_cast<result_type>(Certified_approx_predicate(std::forward(a)...));
}
#endif

} // end namespace CGAL

#endif // CGAL_UNFILTERED_PREDICATE_ADAPTOR_H

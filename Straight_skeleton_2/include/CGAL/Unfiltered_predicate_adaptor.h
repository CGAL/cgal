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

  template <class A1>
  result_type
  operator()(const A1 &a1) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    return static_cast<result_type>(Certified_approx_predicate(a1));
  }
#endif

  template <class A1, class A2>
  result_type
  operator()(const A1 &a1, const A2 &a2) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    return static_cast<result_type>(Certified_approx_predicate(a1, a2));
  }
#endif

  template <class A1, class A2, class A3>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3));
  }
#endif

  template <class A1, class A2, class A3, class A4>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3,a4));
  }
#endif

  template <class A1, class A2, class A3, class A4, class A5>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3, a4, a5));
  }
#endif

  template <class A1, class A2, class A3, class A4, class A5, class A6>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3, a4, a5, a6));
  }
#endif

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3, a4, a5, a6, a7));
  }
#endif

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3, a4, a5, a6, a7, a8));
  }
#endif

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8,
             const A9 &a9 ) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3, a4, a5, a6, a7, a8, a9));
  }
#endif

  // Idem for more than 9 arguments.  Do it on demand.
};

#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
template <class CAP>
  template <class A1>
typename Unfiltered_predicate_adaptor<CAP>::result_type
Unfiltered_predicate_adaptor<CAP>::
  operator()(const A1 &a1) const
{
  return static_cast<result_type>(Certified_approx_predicate(a1));
}

template <class CAP>
  template <class A1, class A2>
typename Unfiltered_predicate_adaptor<CAP>::result_type
Unfiltered_predicate_adaptor<CAP>::
  operator()(const A1 &a1, const A2 &a2) const
{
  return static_cast<result_type>(Certified_approx_predicate(a1, a2));
}

template <class CAP>
  template <class A1, class A2, class A3>
typename Unfiltered_predicate_adaptor<CAP>::result_type
Unfiltered_predicate_adaptor<CAP>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
{
  return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3));
}

template <class CAP>
  template <class A1, class A2, class A3, class A4>
typename Unfiltered_predicate_adaptor<CAP>::result_type
Unfiltered_predicate_adaptor<CAP>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
{
  return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3, a4));
}

template <class CAP>
  template <class A1, class A2, class A3, class A4, class A5>
typename Unfiltered_predicate_adaptor<CAP>::result_type
Unfiltered_predicate_adaptor<CAP>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5) const
{
  return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3, a4, a5));
}

template <class CAP>
  template <class A1, class A2, class A3, class A4, class A5, class A6>
typename Unfiltered_predicate_adaptor<CAP>::result_type
Unfiltered_predicate_adaptor<CAP>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6) const
{
  return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3, a4, a5, a6));
}

template <class CAP>
  template <class A1, class A2, class A3, class A4, class A5, class A6, class A7>
typename Unfiltered_predicate_adaptor<CAP>::result_type
Unfiltered_predicate_adaptor<CAP>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7 ) const
{
  return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3, a4, a5, a6, a7));
}

template <class CAP>
  template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8>
typename Unfiltered_predicate_adaptor<CAP>::result_type
Unfiltered_predicate_adaptor<CAP>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8 ) const
{
  return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3, a4, a5, a6, a7, a8));
}

template <class CAP>
  template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9>
typename Unfiltered_predicate_adaptor<CAP>::result_type
Unfiltered_predicate_adaptor<CAP>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8, 
             const A9 &a9 ) const
{
  return static_cast<result_type>(Certified_approx_predicate(a1, a2, a3, a4, a5, a6, a7, a8, a9));
}
#endif

} // end namespace CGAL

#endif // CGAL_UNFILTERED_PREDICATE_ADAPTOR_H

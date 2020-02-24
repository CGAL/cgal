// Copyright (c) 2007 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// 
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>

#ifndef CGAL_STRAIGHT_SKELETON_ASSERTIONS_H
#define CGAL_STRAIGHT_SKELETON_ASSERTIONS_H

#include <CGAL/license/Straight_skeleton_2.h>


#if defined(CGAL_STRAIGHT_SKELETON_NO_POSTCONDITIONS) \
  || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_STRAIGHT_SKELETON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_stskel_expensive_postcondition(EX)         (static_cast<void>(0))
#  define CGAL_stskel_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_stskel_expensive_postcondition_code(CODE)
#else
#  define CGAL_stskel_expensive_postcondition(EX)         (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_stskel_expensive_postcondition_msg(EX,MSG) (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_stskel_expensive_postcondition_code(CODE)  CODE
#endif 


namespace CGAL {

namespace {

template<class Handle> inline bool handle_assigned ( Handle const& aH )
{
  Handle null ;
  return aH != null ;
}

}

} // end namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_ASSERTIONS_H //
// EOF //
 

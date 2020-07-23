// Copyright (c) 2013
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot <sebastien.loriot@cgal.org>


#ifndef CGAL_RESULT_OF_H
#define CGAL_RESULT_OF_H

#include <CGAL/config.h>
#include <CGAL/disable_warnings.h>
#include <type_traits>

// Address the warning C4003: not enough actual parameters for macro 'BOOST_PP_SEQ_DETAIL_IS_NOT_EMPTY'
// result_of.hpp includes files from boost/preprocessor
// This concerns boost 1_65_1
#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable: 4003)
#endif
#include <boost/utility/result_of.hpp>
#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#include <boost/version.hpp>

namespace CGAL{

namespace cpp11{

// When we switch to C++17, we should use the std::invoke_result instead.
// #if __cplusplus >= 201703L
//   template<typename F>
//   struct result_of : public std::invoke_result<F>
//   {

//   };
// #else
  using std::result_of;
// #endif

}

}

#include <CGAL/enable_warnings.h>

#endif //CGAL_RESULT_OF_H

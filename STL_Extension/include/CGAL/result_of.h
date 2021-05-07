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

#define CGAL_DEPRECATED_HEADER "<CGAL/result_of.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/config.h>"

#include <CGAL/config.h>
#include <CGAL/disable_warnings.h>

#if CGAL_CXX20 || __cpp_lib_is_invocable>=201703L

  // C++>=17

#elif CGAL_CXX11

  #include <type_traits>

#else // C++<11

  // Address the warning C4003: not enough actual parameters for macro 'BOOST_PP_SEQ_DETAIL_IS_NOT_EMPTY'
  // result_of.hpp includes files from boost/preprocessor
  // This concerns boost 1_65_1
  #if defined(BOOST_MSVC)
    #pragma warning(push)
    #pragma warning(disable: 4003)
  #endif
    #include <boost/utility/result_of.hpp>
  #if defined(BOOST_MSVC)
    #pragma warning(pop)
  #endif
  #include <boost/version.hpp>

#endif // end C++<11

namespace CGAL {
namespace cpp11 {

#if CGAL_CXX20 || __cpp_lib_is_invocable>=201703L

  template<typename Signature> class result_of;
  template<typename F, typename... Args>
  class result_of<F(Args...)> : public std::invoke_result<F, Args...> { };

#elif CGAL_CXX11

  using std::result_of;

#else // C++<11

  using boost::result_of;

#endif // end C++<11

} // end cpp11
} // end CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_RESULT_OF_H

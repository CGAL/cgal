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

// Even if for now we use the tr1 result_of implementation, we use the cpp11
// namespace since in the future, that's the decltype version that will be used
namespace cpp11{

template<typename F>
struct result_of
{
  // from boost 1.44 release note https://www.boost.org/users/history/version_1_44_0.html :
  //    New template boost::tr1_result_of that implements the TR1 ResultOf protocol even if boost::result_of uses the C++0x decltype-based implementation.
  #if BOOST_VERSION < 104400
  typedef typename boost::result_of<F>::type type;
  #else
  typedef typename boost::tr1_result_of<F>::type type;
  #endif
};

}

}

#include <CGAL/enable_warnings.h>

#endif //CGAL_RESULT_OF_H

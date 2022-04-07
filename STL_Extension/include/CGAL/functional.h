// Copyright (c) 2017
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_FUNCTIONAL_H
#define CGAL_FUNCTIONAL_H

namespace CGAL {

namespace cpp98 {

  /// Replacement for `std::unary_function` that is deprecated since C++11,
  /// and removed from C++17.
  template < typename ArgumentType, typename ResultType>
  struct unary_function {
    typedef ArgumentType argument_type;
    typedef ResultType result_type;
  };

  /// Replacement for `std::binary_function` that is deprecated since C++11,
  /// and removed from C++17.
  template < typename Arg1, typename Arg2, typename Result>
  struct binary_function {
    typedef Arg1 first_argument_type;
    typedef Arg2 second_argument_type;
    typedef Result result_type;
  };

} // namespace cpp98

} // namespace CGAL

#endif // CGAL_FUNCTIONAL_H

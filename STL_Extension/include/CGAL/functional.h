// Copyright (c) 2017  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_UNARY_FUNCTION_H
#define CGAL_UNARY_FUNCTION_H

namespace CGAL {

  /// Replacement for `std::unary_function` that is deprecated since C++11,
  /// and removed from C++17
  template < typename ArgumentType, typename ResultType>
  struct unary_function {
    typedef ArgumentType argument_type;
    typedef ResultType result_type;
  };


  /// Replacement for `std::binary_function` that is deprecated since C++11,
  /// and removed from C++17
  template < typename Arg1, typename Arg2, typename Result>
  struct binary_function {
    typedef Arg1 first_argument_type;
    typedef Arg2 second_argument_type;
    typedef Result result_type;
  };
  
} // namespace CGAL

#endif // CGAL_UNARY_FUNCTION_H

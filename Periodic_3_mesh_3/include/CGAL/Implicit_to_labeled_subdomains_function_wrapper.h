// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mikhail Bogdanov
//
#ifndef CGAL_PERIODIC_3_MESH_3_IMPLICIT_TO_LABELED_SUBDOMAINS_FUNCTION_WRAPPER_H
#define CGAL_PERIODIC_3_MESH_3_IMPLICIT_TO_LABELED_SUBDOMAINS_FUNCTION_WRAPPER_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4180) // qualifier applied to function type has no meaning; ignored
#endif

namespace CGAL {

template<class Function_, class BGT>
class Implicit_to_labeled_subdomains_function_wrapper
{
public:
  // Types
  typedef int                     return_type;
  typedef typename BGT::Point_3   Point_3;

  /// Constructor
  Implicit_to_labeled_subdomains_function_wrapper(Function_& f)
    : r_f_(f)
  { }

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~Implicit_to_labeled_subdomains_function_wrapper() { }

  /// Operator ()
  return_type operator()(const Point_3& p, const bool = true) const
  {
    return ( (r_f_(p)<0) ? 1 : 2 );
  }

private:
  /// Function to wrap
  Function_& r_f_;

};  // end class Implicit_to_labeled_subdomains_function_wrapper

}  // end namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_PERIODIC_3_MESH_3_IMPLICIT_TO_LABELED_SUBDOMAINS_FUNCTION_WRAPPER_H

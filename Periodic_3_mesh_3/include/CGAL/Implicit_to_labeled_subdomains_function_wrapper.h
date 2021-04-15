// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mikhail Bogdanov
//
#ifndef CGAL_PERIODIC_3_MESH_3_IMPLICIT_TO_LABELED_SUBDOMAINS_FUNCTION_WRAPPER_H
#define CGAL_PERIODIC_3_MESH_3_IMPLICIT_TO_LABELED_SUBDOMAINS_FUNCTION_WRAPPER_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <boost/type_traits/is_function.hpp>
#include <boost/mpl/if.hpp>

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
  Implicit_to_labeled_subdomains_function_wrapper(Function_ f) : f_(f)
  { }

  /// Destructor
  ~Implicit_to_labeled_subdomains_function_wrapper() { }

  /// Operator ()
  return_type operator()(const Point_3& p) const
  {
    // here is the important part : both > 0 ---> both the interior and the exterior are meshed
    return ( (f_(p)<0) ? 1 : 2 );
  }

private:
  typedef typename boost::mpl::if_<boost::is_function<Function_>,
                                   Function_*,
                                   Function_>::type Stored_function;

  /// Function to wrap
  Stored_function f_;
};

}  // end namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_PERIODIC_3_MESH_3_IMPLICIT_TO_LABELED_SUBDOMAINS_FUNCTION_WRAPPER_H

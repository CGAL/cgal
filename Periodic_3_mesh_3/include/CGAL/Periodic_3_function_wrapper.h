// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_PERIODIC_3_MESH_3_PERIODIC_3_FUNCTION_WRAPPER_H
#define CGAL_PERIODIC_3_MESH_3_PERIODIC_3_FUNCTION_WRAPPER_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/internal/canonicalize_helper.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>

#include <boost/type_traits/is_function.hpp>
#include <boost/mpl/if.hpp>

namespace CGAL {

template<class Function_, typename BGT_>
class Periodic_3_function_wrapper
{
public:
  // Types
  typedef typename BGT_::FT                                     FT;
  typedef typename BGT_::Point_3                                Point;
  typedef typename BGT_::Point_3                                Point_3;
  typedef typename BGT_::Iso_cuboid_3                           Iso_cuboid_3;

  typedef FT                                                    return_type;

  typedef typename details::Periodic_3_mesh_geom_traits_generator<BGT_>::type Geom_traits;

  /// Constructor
  Periodic_3_function_wrapper() { }
  Periodic_3_function_wrapper(Function_ f, const Iso_cuboid_3& domain) : f_(f), gt_(domain) { }

  /// Operator ()
  return_type operator()(const Point_3& p) const {
    return f_(P3T3::internal::robust_canonicalize_point(p, gt_));
  }

private:
  typedef typename boost::mpl::if_<boost::is_function<Function_>,
                                   Function_*,
                                   Function_>::type             Stored_function;

  /// Function to wrap
  Stored_function f_;
  Geom_traits gt_;
};

} // namespace CGAL

#endif // CGAL_PERIODIC_3_MESH_3_PERIODIC_3_FUNCTION_WRAPPER_H

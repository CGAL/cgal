// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
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

#ifndef CGAL_BOOST_PARAMETER_H
#define CGAL_BOOST_PARAMETER_H

#ifdef BOOST_PARAMETER_MAX_ARITY
#  if (BOOST_PARAMETER_MAX_ARITY < 12)
#    error "BOOST_PARAMETER_MAX_ARITY must be at least 12 for CGAL::Mesh_3"
#  endif
#else
#  define  BOOST_PARAMETER_MAX_ARITY 12
#endif
#include <boost/parameter.hpp>
#include <boost/parameter/name.hpp>

namespace CGAL
{
namespace parameters
{
  BOOST_PARAMETER_NAME( (time_limit, tag) time_limit_ )
  BOOST_PARAMETER_NAME( (convergence, tag) convergence_)
  BOOST_PARAMETER_NAME( (max_iteration_number, tag) max_iteration_number_ )
  BOOST_PARAMETER_NAME( (freeze_bound, tag) freeze_bound_)
} // parameters
} // CGAL

#endif // CGAL_BOOST_PARAMETER_H

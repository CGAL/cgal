// Copyright (c) 2009-2018 GeometryFactory (France).
// All rights reserved.
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
// Author(s) : Jane Tournois, Laurent Rineau
//

#ifndef CGAL_MESHES_GLOBAL_PARAMETERS_H
#define CGAL_MESHES_GLOBAL_PARAMETERS_H

#include <CGAL/disable_warnings.h>

#include <CGAL/config.h>

#include <boost/parameter.hpp>
#include <boost/parameter/name.hpp>

#if ( defined( __clang__ ) || (BOOST_GCC >= 40600 ) ) && (BOOST_VERSION < 106000)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

namespace CGAL
{
namespace parameters
{

BOOST_PARAMETER_NAME( (max_iteration_number, tag) max_iteration_number_ )
BOOST_PARAMETER_NAME( (convergence, tag) convergence_)
BOOST_PARAMETER_NAME( (time_limit, tag) time_limit_ )
BOOST_PARAMETER_NAME( (freeze_bound, tag) freeze_bound_)

}//end namespace parameters
}//end namespace CGAL

//CGAL_PRAGMA_DIAG_POP
#if ( defined( __clang__ ) || (BOOST_GCC >= 40600 ) ) && (BOOST_VERSION < 106000)
#pragma GCC diagnostic pop
#endif


#include <CGAL/enable_warnings.h>

#endif

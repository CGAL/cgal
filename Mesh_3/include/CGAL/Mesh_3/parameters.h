// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_MESH_3_PARAMETERS_H
#define CGAL_MESH_3_PARAMETERS_H

#include <CGAL/license/Mesh_3.h>
#include <CGAL/Mesh_error_code.h>
#include <CGAL/Mesh_3/sliver_criteria.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/STL_Extension/internal/Has_features.h>

namespace CGAL {
namespace parameters {

#define CGAL_NP_BASE internal_np::No_property
#define CGAL_NP_BUILD(P, V) P(V)

#include <CGAL/STL_Extension/internal/mesh_parameters_interface.h>

#undef CGAL_NP_BASE
#undef CGAL_NP_BUILD

} } // end of CGAL::parameters namespace


#endif //CGAL_MESH_3_PARAMETERS_H

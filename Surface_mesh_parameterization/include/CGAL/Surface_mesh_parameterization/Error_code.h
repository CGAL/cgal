// Copyright (c) 2016  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_ERROR_CODE_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_ERROR_CODE_H

#include <CGAL/license/Surface_mesh_parameterization.h>

/// \file Error_code.h

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup PkgSurfaceMeshParameterizationEnums
///
/// List of errors detected by this package
enum Error_code
{
  OK,                             ///< Success
  ERROR_EMPTY_MESH,               ///< Input mesh is empty
  ERROR_NON_TRIANGULAR_MESH,      ///< Input mesh is not triangular
  ERROR_NO_TOPOLOGICAL_DISC,      ///< Input mesh is not a topological disc
  ERROR_NO_TOPOLOGICAL_BALL,      ///< Input mesh is not a topological ball
  ERROR_BORDER_TOO_SHORT,         ///< This border parameterization requires a longer border
  ERROR_NON_CONVEX_BORDER,        ///< This parameterization method requires a convex border
  ERROR_CANNOT_SOLVE_LINEAR_SYSTEM,///< Cannot solve linear system
  ERROR_NO_1_TO_1_MAPPING,        ///< Parameterization failed: no one-to-one mapping
  ERROR_WRONG_PARAMETER           ///< A method received an unexpected parameter
};

/// \ingroup PkgSurfaceMeshParameterizationEnums
/// \brief gets the message corresponding to an error code.
/// \param error_code The code returned by `parameterize()`
/// \return The string describing the error code.
inline const char* get_error_message(int error_code)
{
  // Messages corresponding to Error_code list above. Must be kept in sync!
  static const char* error_message[ERROR_WRONG_PARAMETER+1] = {
    "Success",
    "Input mesh is empty",
    "Input mesh is not triangular",
    "Input mesh is not a topological disc",
    "This border parameterization requires a longer border",
    "This parameterization method requires a convex border",
    "Cannot solve linear system",
    "Parameterization failed: no one-to-one mapping",
    "A method received an unexpected parameter"
  };

  if(error_code > ERROR_WRONG_PARAMETER || error_code < 0)
    return "Unknown error";
  else
    return error_message[error_code];
}

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_ERROR_CODE_H

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

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_ORBIFOLD_ENUMS_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_ORBIFOLD_ENUMS_H

#include <CGAL/license/Surface_mesh_parameterization.h>

/// \file orbifold_enums.h

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup PkgSurfaceMeshParameterizationEnums
///
/// The two possible weight types available in the Orbifold Tutte parameterization.
enum Weight_type
{
  Cotangent = 0, ///< When Cotangent weights are used, the orbifold-Tutte embedding
                 /// globally minimizes the Dirichlet energy and approximates conformal mappings.
  Mean_value     ///< Mean Value Coordinate weights are guaranteed to generate positive edge weights,
                 /// and the parameterization is guaranteed to be injective.
};

/// \ingroup PkgSurfaceMeshParameterizationEnums
///
/// A classification type for the cones used in Orbifold Tutte parameterization.
enum Cone_type
{
  First_unique_cone = 0, ///< Marker for the cone found at the beginning of the seam.
  Second_unique_cone,    ///< Marker for the cone found at the end of the seam.
  Duplicated_cone        ///< Marker for all the other cones. Cones are duplicated in the sense
                         /// that when the seam is "opened", the cone appears
                         /// at two different positions.
};

/// \ingroup PkgSurfaceMeshParameterizationEnums
///
/// The four Orbifold types available in the Orbifold Tutte parameterization.
/// The different shapes result from the number of cones and the angle constraints
/// at the cones.
enum Orbifold_type
{
  Square = 0,   ///< Three cones, forming a square-shaped basic tile.
  Diamond,      ///< Three cones, forming a diamond-shaped basic tile.
  Triangle,     ///< Three cones, forming a triangle-shaped basic tile.
  Parallelogram ///< Four cones, forming a parallelogram-shaped basic tile.
};

/// \ingroup PkgSurfaceMeshParameterizationEnums
/// \brief Convert the orbifold type to a literal message.
/// \param orb_type the integer value in the enum
/// \return the string describing the Orbifold type.
const char* get_orbifold_type(int orb_type)
{
  // Messages corresponding to the different orbifold types.
  static const char* type[Parallelogram+1] = {
    "Square",
    "Diamond",
    "Triangle",
    "Parallelogram"
  };

  if(orb_type > Parallelogram || orb_type < 0)
    return "Unknown orbifold type";
  else
    return type[orb_type];
}

} // namespace Surface_mesh_parameterization
} // namespace CGAL

#endif

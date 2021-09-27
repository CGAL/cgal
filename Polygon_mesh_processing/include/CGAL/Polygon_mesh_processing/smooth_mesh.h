// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé
//                 Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTH_MESH_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTH_MESH_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Polygon_mesh_processing/smooth_mesh.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>"
#include <CGAL/Installation/internal/deprecation_warning.h>

#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>

#ifndef CGAL_NO_DEPRECATED_CODE

#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
* \ingroup PMPDeprecated
*
* \deprecated This function is deprecated since \cgal 5.5,
* `CGAL::angle_and_area_smoothing()` should be used instead.
*/
template<typename TriangleMesh, typename FaceRange, typename NamedParameters>
CGAL_DEPRECATED void smooth_mesh(const FaceRange& faces,
                 TriangleMesh& tmesh,
                 const NamedParameters& np)
{
  angle_and_area_smoothing(faces, tmesh, np);
}

///\cond SKIP_IN_MANUAL
template <typename FaceRange, typename TriangleMesh>
CGAL_DEPRECATED void smooth_mesh(const FaceRange& face_range, TriangleMesh& tmesh)
{
  smooth_mesh(face_range, tmesh, parameters::all_default());
}

template <typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED void smooth_mesh(TriangleMesh& tmesh, const CGAL_PMP_NP_CLASS& np)
{
  smooth_mesh(faces(tmesh), tmesh, np);
}

template<typename TriangleMesh>
CGAL_DEPRECATED void smooth_mesh(TriangleMesh& tmesh)
{
  smooth_mesh(faces(tmesh), tmesh, parameters::all_default());
}
///\endcond SKIP_IN_MANUAL

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif //#ifndef CGAL_NO_DEPRECATED_CODE

#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTH_MESH_H

// Copyright (c) 2025
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dema Nasrawt and Aidan Pawlak

#ifndef CGAL_IO_READ_GLTF_H
#define CGAL_IO_READ_GLTF_H

#include <CGAL/IO/GLTF/read_gltf.h>
#include <string>

namespace CGAL {
namespace IO {

/// \ingroup PkgStreamSupportIoFuncsGLTF
///
/// \brief reads the content of the file `fname` into `points` and `polygons`, using the \ref IOStreamGLTF.
///
/// \attention The polygon soup is not cleared, and the data from the file are appended.
///
/// \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
/// \tparam PolygonRange a model of the concepts `SequenceContainer` and `BackInsertionSequence`
///                      whose `value_type` is itself a model of the concept `SequenceContainer`
///                      and `BackInsertionSequence` whose `value_type` is an unsigned integer type
///                      convertible to `std::size_t`
///
/// \param filename the path to the input file
/// \param points points of the soup of polygons
/// \param polygons a range of polygons. Each element in it describes a polygon
///        using the indices of the points in `points`.
///
/// \returns `true` if the reading was successful, `false` otherwise.
template <typename PointRange, typename PolygonRange>
bool read_GLTF(const std::string& filename,
               PointRange& points,
               PolygonRange& polygons) {
    return internal::read_GLTF(filename, points, polygons);
}

// Public entry point that wraps your internal implementation
template <typename PolygonMesh>
bool read_GLTF(const std::string& filename, PolygonMesh& mesh) {
    return internal::read_GLTF(filename, mesh);
}

} // namespace IO
} // namespace CGAL

#endif

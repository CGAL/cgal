// Copyright (c) 2015  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Lutz Kettner
//             Andreas Fabri
//             Maxime Gimeno

#ifndef CGAL_IO_OBJ_H
#define CGAL_IO_OBJ_H

#include <CGAL/IO/OBJ/File_writer_wavefront.h>
#include <CGAL/IO/OBJ/OBJ_reader.h>

namespace CGAL {

//! \ingroup IOstreamFunctions
//!
///reads a file in .obj format.
///
///\tparam Points_3 a RandomAccessContainer of Point_3,
///\tparam Faces a RandomAccessContainer of RandomAccessContainer of std::size_t
///
/// \see IOStreamOBJ
template <class Points_3, class Faces>
bool read_OBJ(std::istream& input, Points_3 &points, Faces &faces);

} // namespace CGAL

#endif // CGAL_IO_OBJ_H

// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)    : Samuel Hornus

#ifndef CGAL_INTERNAL_TRIANGULATION_DUMMY_TDS_H
#define CGAL_INTERNAL_TRIANGULATION_DUMMY_TDS_H

#include <CGAL/license/Triangulation.h>


namespace CGAL {

namespace internal {
namespace Triangulation {

struct Dummy_TDS
{ 
    struct Vertex {};
    struct Vertex_handle {};
    struct Vertex_iterator {};
    struct Vertex_const_handle {};
    struct Vertex_const_iterator {};
    struct Full_cell {};
    struct Full_cell_handle {};
    struct Full_cell_iterator {};
    struct Full_cell_const_handle {};
    struct Full_cell_const_iterator {};
    struct Vertex_handle_const_iterator {};
    struct Full_cell_data {};
};

} // namespace Triangulation
} // namespace internal

} //namespace CGAL

#endif // CGAL_INTERNAL_TRIANGULATION_DUMMY_TDS_H

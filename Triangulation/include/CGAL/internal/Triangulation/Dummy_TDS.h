// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

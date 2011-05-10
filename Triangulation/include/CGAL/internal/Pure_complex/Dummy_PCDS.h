// Copyright (c) 2009 INRIA Sophia-Antipolis (France),
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)    : Samuel Hornus

#ifndef CGAL_INTERNAL_TRIANGULATION_DUMMY_TDS_H
#define CGAL_INTERNAL_TRIANGULATION_DUMMY_TDS_H

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
    struct Simplex {};
    struct Simplex_handle {};
    struct Simplex_iterator {};
    struct Simplex_const_handle {};
    struct Simplex_const_iterator {};
    struct Vertex_handle_const_iterator {};
};

}; // namespace Triangulation
}; // namespace internal

} //namespace CGAL

#endif // CGAL_INTERNAL_TRIANGULATION_DUMMY_TDS_H

// Copyright (c) 1997  ETH Zurich (Switzerland).
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
// 
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>)

#ifndef CGAL_POLYHEDRON_TRAITS_WITH_NORMALS_3_H
#define CGAL_POLYHEDRON_TRAITS_WITH_NORMALS_3_H 1

#include <CGAL/basic.h>

namespace CGAL {

template < class Kernel_ >
class Polyhedron_traits_with_normals_3 {
public:
    typedef Kernel_                   Kernel;
    typedef typename Kernel::Point_3  Point_3;
    typedef typename Kernel::Vector_3 Plane_3;

    typedef typename Kernel::Construct_opposite_vector_3 
                                      Construct_opposite_plane_3;
private:
    Kernel m_kernel;

public:
    Polyhedron_traits_with_normals_3() {}
    Polyhedron_traits_with_normals_3( const Kernel& kernel)
        : m_kernel(kernel) {}

    Construct_opposite_plane_3 construct_opposite_plane_3_object() const {
        return m_kernel.construct_opposite_vector_3_object();
    }
};

} //namespace CGAL

#endif // CGAL_POLYHEDRON_TRAITS_WITH_NORMALS_3_H //
// EOF //

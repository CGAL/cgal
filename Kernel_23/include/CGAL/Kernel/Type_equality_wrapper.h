// Copyright (c) 2003  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_TYPE_EQUALITY_WRAPPER_H
#define CGAL_KERNEL_TYPE_EQUALITY_WRAPPER_H

#include <CGAL/user_classes.h>

namespace CGAL {

// This is a kernel wrapper which provides the type equality between
// Kernel::Point_2 and CGAL::Point_2<Kernel>, by deriving from
// K_base::Point_2 (and similar for the other types).

template < typename K_base, typename Kernel_ >
struct Type_equality_wrapper
  : public K_base
{
    typedef K_base                                  Kernel_base;

#define CGAL_Kernel_obj(X)   typedef CGAL::X<Kernel_> X;

#include <CGAL/Kernel/interface_macros.h>

    // Undocumented stuff.
    typedef CGAL::Conic_2<Kernel_>                   Conic_2;
    typedef CGAL::Aff_transformation_2<Kernel_>      Aff_transformation_2;
    typedef CGAL::Aff_transformation_3<Kernel_>      Aff_transformation_3;
};

} //namespace CGAL

#endif // CGAL_KERNEL_TYPE_EQUALITY_WRAPPER_H

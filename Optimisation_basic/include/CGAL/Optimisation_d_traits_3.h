// Copyright (c) 1997-2001  
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
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_OPTIMISATION_D_TRAITS_3_H
#define CGAL_OPTIMISATION_D_TRAITS_3_H

// includes
#  include <CGAL/Optimisation/Access_dimension_3.h>
#  include <CGAL/Optimisation/Access_coordinates_begin_3.h>
#  include <CGAL/Optimisation/Construct_point_3.h>

namespace CGAL {

// Class declaration
// =================
template < class K_, class ET_ = typename K_::RT,
                     class NT_ = typename K_::RT >
class Optimisation_d_traits_3;

// Class interface
// ===============
template < class K_, class ET_, class NT_>
class Optimisation_d_traits_3 {
  public:
    // self
    typedef  K_                         K;
    typedef  ET_                        ET;
    typedef  NT_                        NT;
    typedef  Optimisation_d_traits_3<K,ET,NT>
                                        Self;

    // types
    typedef  typename K::Point_3        Point_d;

    typedef  typename K::Rep_tag        Rep_tag;

    typedef  typename K::RT             RT;
    typedef  typename K::FT             FT;

    typedef  CGAL::Access_dimension_3<K>      Access_dimension_d;
    typedef  CGAL::Access_coordinates_begin_3<K>
                                        Access_coordinates_begin_d;

    typedef  _Construct_point_3<K>       Construct_point_d;
    // this does not (yet) work:
    // typedef  typename K::Construct_point_3    Construct_point_d;

    // creation
    Optimisation_d_traits_3( ) { }
    Optimisation_d_traits_3( const Optimisation_d_traits_3<K_,ET_,NT_>&) {}

    // operations
    Access_dimension_d
    access_dimension_d_object( ) const
        { return Access_dimension_d(); }

    Access_coordinates_begin_d
    access_coordinates_begin_d_object( ) const
        { return Access_coordinates_begin_d(); }

    Construct_point_d
    construct_point_d_object( ) const
        { return Construct_point_d(); }
};

} //namespace CGAL

#endif // CGAL_OPTIMISATION_D_TRAITS_3_H

// ===== EOF ==================================================================

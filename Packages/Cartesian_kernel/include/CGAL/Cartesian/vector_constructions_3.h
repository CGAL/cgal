// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_VECTOR_CONSTRUCTIONS_3_H
#define CGAL_CARTESIAN_VECTOR_CONSTRUCTIONS_3_H

#include <CGAL/Cartesian/Vector_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class K >
inline
typename K::Vector_3
cross_product(const VectorC3<K> &v, const VectorC3<K> &w)
{
  return K().construct_cross_product_vector_3_object()(v, w);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_VECTOR_CONSTRUCTIONS_3_H

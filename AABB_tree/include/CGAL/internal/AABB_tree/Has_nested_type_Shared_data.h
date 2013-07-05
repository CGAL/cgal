// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://sloriot@scm.gforge.inria.fr/svn/cgal/branches/features/AABB_tree-one_primitive_per_object-sloriot/AABB_tree/include/CGAL/AABB_triangle_primitive.h $
// $Id: AABB_triangle_primitive.h 69127 2012-05-14 16:10:00Z sloriot $
//
//
// Author(s)     : Sebastien Loriot
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_INTERNAL_AABB_TREE_HAS_NESTED_TYPE_SHARED_DATA_H
#define CGAL_INTERNAL_AABB_TREE_HAS_NESTED_TYPE_SHARED_DATA_H

#include <boost/mpl/has_xxx.hpp>

namespace CGAL{

namespace internal{

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Shared_data,Shared_data,false)
  
} } //namespace CGAL

#endif //CGAL_INTERNAL_AABB_TREE_HAS_NESTED_TYPE_SHARED_DATA_H

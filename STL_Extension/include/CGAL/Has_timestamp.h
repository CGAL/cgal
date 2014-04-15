// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// Copyright (c) 2014 GeometryFactory Sarl (France)
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
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_HAS_TIMESTAMP_H
#define CGAL_HAS_TIMESTAMP_H

#include <boost/mpl/has_xxx.hpp>
#include <CGAL/tags.h>

namespace CGAL {

namespace internal {

  BOOST_MPL_HAS_XXX_TRAIT_DEF(Has_timestamp)

  // to have Mesh_3 deterministic,
  // a partial specialization of this class should be written next to
  // every class that implements concepts MeshCellBase_3 or MeshVertexBase_3
  template <typename T, bool has_ts = has_Has_timestamp<T>::value>
  struct Has_timestamp : public T::Has_timestamp
    // when T does not have a partial specialization of Has_timestamp
  {};

  template <typename T>
  struct Has_timestamp<T, false> : public Tag_false
    // when T does not have a partial specialization of Has_timestamp
  {};

} // end namespace internal
} // end namespace CGAL

#endif // CGAL_HAS_TIMESTAMP_H

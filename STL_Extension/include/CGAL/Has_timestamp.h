// Copyright (c) 2014 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jane Tournois

#ifndef CGAL_HAS_TIMESTAMP_H
#define CGAL_HAS_TIMESTAMP_H

#include <boost/mpl/has_xxx.hpp>
#include <CGAL/tags.h>

namespace CGAL {

namespace internal {

  BOOST_MPL_HAS_XXX_TRAIT_DEF(Has_timestamp)

  // Used by Compact container to make the comparison of iterator
  // depending on the insertion order rather than the object address
  // when the object class defines a Has_timestamp tag
  // This is for example used in to make Mesh_3 deterministic, see
  // classes implementing concepts MeshCellBase_3 and MeshVertexBase_3
  template <typename T, bool has_ts = has_Has_timestamp<T>::value>
  struct Has_timestamp : public T::Has_timestamp
    // when T has a Has_timestamp tag
  {};

  template <typename T>
  struct Has_timestamp<T, false> : public Tag_false
    // when T does not have a Has_timestamp tag
  {};

} // end namespace internal
} // end namespace CGAL

#endif // CGAL_HAS_TIMESTAMP_H

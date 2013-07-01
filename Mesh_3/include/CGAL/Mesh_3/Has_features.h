// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// Copyright (c) 2011 GeometryFactory Sarl (France)
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
// $URL: https://scm.gforge.inria.fr/svn/cgal/branches/features/Mesh_3-experimental-GF/Mesh_3/include/CGAL/Mesh_3/Has_features.h $
// $Id: Has_features.h 61189 2011-02-11 18:09:01Z lrineau $
//
//
// Author(s)     : Stéphane Tayeb, Laurent Rineau

#ifndef CGAL_MESH_3_HAS_FEATURES_H
#define CGAL_MESH_3_HAS_FEATURES_H

#include <boost/mpl/has_xxx.hpp>
#include <CGAL/tags.h>

namespace CGAL {

namespace internal {
namespace Mesh_3 {
  // A type has_Has_features to check if type 'Has_features' is a nested
  // type of any class
  BOOST_MPL_HAS_XXX_TRAIT_DEF(Has_features)

  template <typename Mesh_domain, 
            bool has_Has_features = has_Has_features<Mesh_domain>::value>
  struct Has_features : 
    public CGAL::Boolean_tag<Mesh_domain::Has_features::value>
    // when Mesh_domain has the nested type Has_features
  {};

  template <typename Mesh_domain>
  struct Has_features<Mesh_domain, false> : public CGAL::Tag_false
    // when Mesh_domain does not have the nested type Has_features
  {};
} // end namespace internal::Mesh_3
} // end namespace internal
} // end namespace CGAL

#endif // CGAL_MESH_3_HAS_FEATURES_H
